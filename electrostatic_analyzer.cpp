#include "electrostatic_analyzer.h"

#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <string>
#include <chrono>
#include <vector>

ElectrostaticAnalyzer::ElectrostaticAnalyzer()
	: mesh_node_count(0), mesh_element_count(0), dt(0.0), nstep(200), unit(1.0)
{
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;
	last_solver_iterations = 0;
	last_solver_time = 0.0;
	last_solver_method = "ICCG";
	step_output_prefix = "tempa.dat";
	report_output_path = "analysis_summary.txt";
}

void ElectrostaticAnalyzer::set_output_paths(const std::string &step_prefix, const std::string &report_path)
{
	step_output_prefix = step_prefix.empty() ? std::string("tempa.dat") : step_prefix;
	report_output_path = report_path.empty() ? std::string("analysis_summary.txt") : report_path;
}

void ElectrostaticAnalyzer::add_fixed_node(int node, double value, bool overwrite)
{
	if (node < 0 || node >= (int)fixed_node.size()) return;
	if (!fixed_node[node] || overwrite) {
		fixed_node[node] = true;
		fixed_value[node] = value;
	}
}

void ElectrostaticAnalyzer::build_sparse_from_rows(const std::vector<std::vector<std::pair<int, double>>> &rows,
	std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &values) const
{
	const int n = (int)nodes.size();
	row_ptr.assign(n + 1, 0);
	col_idx.clear();
	values.clear();

	for (int r = 0; r < n; ++r) {
		std::vector<std::pair<int, double>> row = rows[r];
		std::sort(row.begin(), row.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
			return a.first < b.first;
		});

		for (size_t i = 0; i < row.size(); ) {
			int c = row[i].first;
			double sum = 0.0;
			size_t j = i;
			while (j < row.size() && row[j].first == c) {
				sum += row[j].second;
				++j;
			}
			if (std::fabs(sum) > 1e-30) {
				col_idx.push_back(c);
				values.push_back(sum);
			}
			i = j;
		}
		row_ptr[r + 1] = (int)col_idx.size();
	}
}

void ElectrostaticAnalyzer::build_time_step_matrix(double scale_eps,
	std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &values) const
{
	auto t_start = std::chrono::high_resolution_clock::now();
	const int n = (int)nodes.size();
	row_ptr.assign(n + 1, 0);
	col_idx.clear();
	values.clear();
	if (n == 0) {
		return;
	}

	// global_K/global_E are already sorted and aggregated by column per row.
	// Merge them directly to avoid per-step row copies and sorting.
	for (int r = 0; r < n; ++r) {
		int pK = global_K_row_ptr[r];
		int pKEnd = global_K_row_ptr[r + 1];
		int pE = global_E_row_ptr[r];
		int pEEnd = global_E_row_ptr[r + 1];

		while (pK < pKEnd || pE < pEEnd) {
			int c;
			double v;
			if (pE >= pEEnd || (pK < pKEnd && global_K_col_idx[pK] < global_E_col_idx[pE])) {
				c = global_K_col_idx[pK];
				v = global_K_values[pK];
				++pK;
			} else if (pK >= pKEnd || global_E_col_idx[pE] < global_K_col_idx[pK]) {
				c = global_E_col_idx[pE];
				v = global_E_values[pE] * scale_eps;
				++pE;
			} else {
				c = global_K_col_idx[pK];
				v = global_K_values[pK] + global_E_values[pE] * scale_eps;
				++pK;
				++pE;
			}

			if (std::fabs(v) > 1e-30) {
				col_idx.push_back(c);
				values.push_back(v);
			}
		}
		row_ptr[r + 1] = (int)col_idx.size();
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	printf("build_time_step_matrix (scale=%.6e): %.6f s\n", scale_eps, elapsed);
}

void ElectrostaticAnalyzer::build_sparse_matrix()
{
	auto t_start = std::chrono::high_resolution_clock::now();
	build_sparse_from_rows(global_K_rows, global_K_row_ptr, global_K_col_idx, global_K_values);
	build_sparse_from_rows(global_E_rows, global_E_row_ptr, global_E_col_idx, global_E_values);

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	printf("build_sparse_matrix: %.6f s\n", elapsed);
}

void ElectrostaticAnalyzer::apply_dirichlet_constraints(std::vector<double> &rhs)
{
	const int n = (int)nodes.size();
	// Single pass over CSR matrix: adjust rhs for fixed columns and zero non-diagonal
	// entries that reference fixed nodes. This is O(nnz).
	for (int r = 0; r < n; ++r) {
		for (int p = global_K_row_ptr[r]; p < global_K_row_ptr[r + 1]; ++p) {
			int c = global_K_col_idx[p];
			if (fixed_node[c]) {
				if (r != c) {
					rhs[r] -= global_K_values[p] * fixed_value[c];
					global_K_values[p] = 0.0;
				}
			}
		}
	}

	// Ensure diagonal entries for fixed nodes are 1.0 and other row entries are zero
	for (int i = 0; i < n; ++i) {
		if (!fixed_node[i]) continue;
		bool has_diag = false;
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			if (global_K_col_idx[p] == i) {
				global_K_values[p] = 1.0;
				has_diag = true;
			} else {
				global_K_values[p] = 0.0;
			}
		}
		if (!has_diag) {
			int insert_pos = global_K_row_ptr[i + 1];
			global_K_col_idx.insert(global_K_col_idx.begin() + insert_pos, i);
			global_K_values.insert(global_K_values.begin() + insert_pos, 1.0);
			for (int r = i + 1; r <= n; ++r) ++global_K_row_ptr[r];
		}
		rhs[i] = fixed_value[i];
	}
}

void ElectrostaticAnalyzer::matvec(const std::vector<double> &x, std::vector<double> &y) const
{
	matvec_csr(global_K_row_ptr, global_K_col_idx, global_K_values, x, y);
}

void ElectrostaticAnalyzer::matvec_csr(const std::vector<int> &row_ptr, const std::vector<int> &col_idx,
	const std::vector<double> &values, const std::vector<double> &x, std::vector<double> &y) const
{
	const int n = (int)nodes.size();
	y.assign(n, 0.0);
	for (int r = 0; r < n; ++r) {
		double sum = 0.0;
		for (int p = row_ptr[r]; p < row_ptr[r + 1]; ++p) {
			sum += values[p] * x[col_idx[p]];
		}
		y[r] = sum;
	}
}

void ElectrostaticAnalyzer::build_incomplete_cholesky()
{
	const int n = (int)nodes.size();
	ic_scale.assign(n, 1.0);
	preconditioner_diag.assign(n, 0.0);
	L_row_ptr.assign(n + 1, 0);
	L_col_idx.clear();
	L_values.clear();

	// Symmetric diagonal scaling to reduce coefficient spread before IC(0).
	const double shift = 1e-30;
	for (int i = 0; i < n; ++i) {
		double diag = 0.0;
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			if (global_K_col_idx[p] == i) {
				diag = global_K_values[p];
				break;
			}
		}
		if (!std::isfinite(diag) || diag <= shift) diag = 1.0;
		preconditioner_diag[i] = diag;
		ic_scale[i] = 1.0 / std::sqrt(diag);
	}

	// build pattern: keep columns <= row (lower triangle including diag)
	for (int i = 0; i < n; ++i) {
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			int c = global_K_col_idx[p];
			if (c <= i) {
				L_col_idx.push_back(c);
				// placeholder for value
				L_values.push_back(0.0);
			}
		}
		L_row_ptr[i + 1] = (int)L_col_idx.size();
	}

	// helper to get A(i,j)
	auto get_A = [&](int row, int col)->double{
		for (int p = global_K_row_ptr[row]; p < global_K_row_ptr[row+1]; ++p) {
			if (global_K_col_idx[p] == col) return ic_scale[row] * global_K_values[p] * ic_scale[col];
		}
		return 0.0;
	};

	// numeric IC(0)
	for (int i = 0; i < n; ++i) {
		// process off-diagonals j < i
		for (int p = L_row_ptr[i]; p < L_row_ptr[i+1]; ++p) {
			int j = L_col_idx[p];
			if (j == i) continue;
			double val = get_A(i, j);

			// compute sum_k L(i,k)*L(j,k) for k < j
			int pi = L_row_ptr[i];
			int pj = L_row_ptr[j];
			while (pi < L_row_ptr[i+1] && pj < L_row_ptr[j+1]) {
				int ki = L_col_idx[pi];
				int kj = L_col_idx[pj];
				if (ki < kj) ++pi;
				else if (ki > kj) ++pj;
				else {
					if (ki < j) {
						val -= L_values[pi] * L_values[pj];
					}
					++pi; ++pj;
				}
			}

			// find L_jj (diagonal of j)
			double Ljj = 0.0;
			for (int q = L_row_ptr[j]; q < L_row_ptr[j+1]; ++q) {
				if (L_col_idx[q] == j) { Ljj = L_values[q]; break; }
			}
			if (std::fabs(Ljj) < 1e-30) Ljj = 1e-30;
			L_values[p] = val / Ljj;
		}

		// diagonal
		double diag = get_A(i, i);
		for (int q = L_row_ptr[i]; q < L_row_ptr[i+1]; ++q) {
			int k = L_col_idx[q];
			if (k == i) continue;
			diag -= L_values[q] * L_values[q];
		}
		if (!std::isfinite(diag) || diag <= 0.0) diag = 1e-30;
		// set diagonal
		for (int q = L_row_ptr[i]; q < L_row_ptr[i+1]; ++q) {
			if (L_col_idx[q] == i) { L_values[q] = std::sqrt(diag); break; }
		}
	}

	// build CSC for L (columns)
	L_csc_col_ptr.assign(n + 1, 0);
	for (size_t idx = 0; idx < L_col_idx.size(); ++idx) {
		int col = L_col_idx[idx];
		++L_csc_col_ptr[col + 1];
	}
	for (int i = 0; i < n; ++i) L_csc_col_ptr[i+1] += L_csc_col_ptr[i];
	L_csc_row_idx.assign(L_col_idx.size(), 0);
	L_csc_values.assign(L_values.size(), 0.0);
	// fill
	std::vector<int> fill_pos = L_csc_col_ptr;
	for (int row = 0; row < n; ++row) {
		for (int p = L_row_ptr[row]; p < L_row_ptr[row+1]; ++p) {
			int col = L_col_idx[p];
			int dst = fill_pos[col]++;
			L_csc_row_idx[dst] = row;
			L_csc_values[dst] = L_values[p];
		}
	}

	ic_built = true;
}

void ElectrostaticAnalyzer::apply_ssor_preconditioner(const std::vector<double> &r, std::vector<double> &z) const
{
	const int n = (int)nodes.size();
	z.assign(n, 0.0);
	std::vector<double> diag_values(n, 1.0);
	const double omega = 1.3;
	const double scale = omega * (2.0 - omega);
	if (preconditioner_diag.size() == (size_t)n) {
		diag_values = preconditioner_diag;
	} else {
		for (int i = 0; i < n; ++i) {
			double diag = 0.0;
			for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
				if (global_K_col_idx[p] == i) {
					diag = global_K_values[p];
					break;
				}
			}
			if (!std::isfinite(diag) || std::fabs(diag) < 1e-30) diag = 1.0;
			diag_values[i] = diag;
		}
	}

	std::vector<double> y(n, 0.0);

	// Forward sweep: (D + L) y = r
	for (int i = 0; i < n; ++i) {
		double diag = diag_values[i];
		if (!std::isfinite(diag) || std::fabs(diag) < 1e-30) diag = 1.0;
		double sum = r[i];
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			int c = global_K_col_idx[p];
			if (c < i) {
				sum -= omega * global_K_values[p] * y[c];
			}
		}
		y[i] = sum / diag;
	}

	// Middle scaling by D
	for (int i = 0; i < n; ++i) {
		y[i] *= diag_values[i];
	}

	// Backward sweep: (D + U) z = D y
	for (int i = n - 1; i >= 0; --i) {
		double diag = diag_values[i];
		if (!std::isfinite(diag) || std::fabs(diag) < 1e-30) diag = 1.0;
		double sum = y[i];
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			int c = global_K_col_idx[p];
			if (c > i) {
				sum -= omega * global_K_values[p] * z[c];
			}
		}
		z[i] = scale * (sum / diag);
	}

	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) z[i] = 0.0;
	}
}

void ElectrostaticAnalyzer::apply_preconditioner(const std::vector<double> &r, std::vector<double> &z) const
{
	const int n = (int)nodes.size();
	z.assign(n, 0.0);
	if (use_ssor_preconditioner) {
		apply_ssor_preconditioner(r, z);
		return;
	}
	if (!ic_built) {
		// fallback: diagonal preconditioner
		for (int i = 0; i < n; ++i) {
			double Aii = 0.0;
			for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i+1]; ++p) if (global_K_col_idx[p] == i) { Aii = global_K_values[p]; break; }
			if (std::fabs(Aii) < 1e-30) z[i] = r[i]; else z[i] = r[i] / Aii;
		}
		return;
	}

	std::vector<double> scaled_r(n, 0.0);
	for (int i = 0; i < n; ++i) {
		scaled_r[i] = r[i] * ic_scale[i];
	}

	// forward solve L * y = r
	std::vector<double> y(n, 0.0);
	for (int i = 0; i < n; ++i) {
		double s = scaled_r[i];
		double Lii = 1.0;
		for (int p = L_row_ptr[i]; p < L_row_ptr[i+1]; ++p) {
			int c = L_col_idx[p];
			if (c == i) { Lii = L_values[p]; break; }
			s -= L_values[p] * y[c];
		}
		if (std::fabs(Lii) < 1e-30) Lii = 1e-30;
		y[i] = s / Lii;
	}

	// backward solve L^T * z = y using CSC
	for (int col = n - 1; col >= 0; --col) {
		double s = y[col];
		double Lii = 1.0;
		// diagonal in CSC is somewhere in column col at row==col
		for (int p = L_csc_col_ptr[col]; p < L_csc_col_ptr[col+1]; ++p) {
			int row = L_csc_row_idx[p];
			if (row == col) { Lii = L_csc_values[p]; break; }
		}
		for (int p = L_csc_col_ptr[col]; p < L_csc_col_ptr[col+1]; ++p) {
			int row = L_csc_row_idx[p];
			if (row == col) continue;
			// L(row,col) contributes to equation for col: subtract L(row,col) * z[row]
			s -= L_csc_values[p] * z[row];
		}
		if (std::fabs(Lii) < 1e-30) Lii = 1e-30;
		z[col] = s / Lii;
	}

	for (int i = 0; i < n; ++i) {
		z[i] *= ic_scale[i];
	}
}

void ElectrostaticAnalyzer::read_mesh(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	int num_nodes = 0, num_elements = 0;
	int dummy1 = 0, dummy2 = 0;
	if (fscanf(fp, "%d %d %d %d %lf", &num_nodes, &num_elements, &dummy1, &dummy2, &unit) != 5) {
		printf("Invalid header in %s\n", filename);
		fclose(fp);
		return;
	}

	elements.assign(num_elements, Element{{0, 0, 0, 0}, 0});
	for (int e = 0; e < num_elements; ++e) {
		int n0 = 0, n1 = 0, n2 = 0, n3 = 0;
		if (fscanf(fp, "%d %d %d %d", &n0, &n1, &n2, &n3) != 4) {
			printf("Invalid element connectivity at %d in %s\n", e + 1, filename);
			fclose(fp);
			return;
		}
		// in.dat は 1-based を想定
		elements[e].nodes[0] = n0 - 1;
		elements[e].nodes[1] = n1 - 1;
		elements[e].nodes[2] = n2 - 1;
		elements[e].nodes[3] = n3 - 1;
	}

	nodes.assign(num_nodes, Node{0.0, 0.0, 0.0});
	for (int i = 0; i < num_nodes; ++i) {
		if (fscanf(fp, "%lf %lf %lf", &nodes[i].x, &nodes[i].y, &nodes[i].z) != 3) {
			printf("Invalid node coordinate at %d in %s\n", i + 1, filename);
			fclose(fp);
			return;
		}
	}

	fclose(fp);
	// store counts for external getters
	mesh_node_count = num_nodes;
	mesh_element_count = num_elements;
	printf("Read mesh: %d nodes, %d elements\n", num_nodes, num_elements);
}

void ElectrostaticAnalyzer::read_material_and_bc(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	int ntblk = 0, ntmblk = 0, ntmate = 0, npt2 = 0, nlv = 0;
	if (fscanf(fp, "%d %d %d %d %d", &ntblk, &ntmblk, &ntmate, &npt2, &nlv) != 5) {
		printf("Invalid header in %s\n", filename);
		fclose(fp);
		return;
	}

	std::vector<std::string> lines;
	char buf[1024];
	while (fgets(buf, sizeof(buf), fp)) {
		std::string line(buf);
		// 空白のみの行を除外
		bool all_space = true;
		for (char c : line) {
			if (!isspace((unsigned char)c)) {
				all_space = false;
				break;
			}
		}
		if (!all_space) lines.push_back(line);
	}
	fclose(fp);

	auto count_tokens = [](const std::string &line) -> int {
		std::istringstream iss(line);
		int count = 0;
		double dummy;
		while (iss >> dummy) ++count;
		return count;
	};

	auto parse_ints = [](const std::string &line, std::vector<int> &out) -> bool {
		std::istringstream iss(line);
		int v;
		out.clear();
		while (iss >> v) out.push_back(v);
		return !out.empty();
	};

	auto parse_doubles = [](const std::string &line, std::vector<double> &out) -> bool {
		std::istringstream iss(line);
		double v;
		out.clear();
		while (iss >> v) out.push_back(v);
		return !out.empty();
	};

	// 先頭行の列数から「境界あり/なし」を推定
	bool has_boundary_section = false;
	for (const std::string &line : lines) {
		int tok = count_tokens(line);
		if (tok == 4) {
			has_boundary_section = true;
		}
		break;
	}

	boundaries.clear();
	material_ranges.clear();
	fixed_materials.clear();
	int idx = 0;
	if (has_boundary_section && ntblk > 0) {
		boundaries.assign(ntblk, BoundaryCondition{0, 0, 0, 0.0});
		for (int i = 0; i < ntblk && idx < (int)lines.size(); ++i, ++idx) {
			std::vector<int> vals;
			std::istringstream iss(lines[idx]);
			double dv = 0.0;
			int iv = 0;
			if (!(iss >> iv)) break;
			int n1 = iv;
			if (!(iss >> iv)) break;
			int n2 = iv;
			if (!(iss >> iv)) break;
			int type = iv;
			if (!(iss >> dv)) break;
			boundaries[i].node1 = n1 - 1;
			boundaries[i].node2 = n2 - 1;
			boundaries[i].type = type;
			boundaries[i].value = dv;
		}
	} else {
		// 境界条件なし
		ntblk = 0;
	}

	// 要素の材料番号割り当て
	for (int i = 0; i < ntmblk && idx < (int)lines.size(); ++i, ++idx) {
		std::vector<int> vals;
		if (!parse_ints(lines[idx], vals) || vals.size() < 3) {
			printf("Invalid material block record at %d in %s\n", i + 1, filename);
			return;
		}
		int start = vals[0];
		int end = vals[1];
		int mat = vals[2];
		int mat_id = mat - 1; // 1-based -> 0-based
		material_ranges.push_back(MaterialBlockRange{start, end, mat_id});
		for (int e = start - 1; e <= end - 1; ++e) {
			if (0 <= e && e < (int)elements.size()) {
				elements[e].material = mat_id;
			}
		}
	}

	// initialize materials with defaults (sigma = 1, epsilon = 1)
	materials.assign(ntmate, Material{1.0,1.0,1.0, 1.0,1.0,1.0, 0.0, 0.0});
	for (int i = 0; i < ntmate && idx < (int)lines.size(); ++i, ++idx) {
		std::vector<double> vals;
		if (!parse_doubles(lines[idx], vals) || (vals.size() != 6 && vals.size() != 4 && vals.size() != 5)) {
			printf("Invalid material record at %d in %s\n", i + 1, filename);
			return;
		}
		if (vals.size() == 6) {
			materials[i].sx = vals[0];
			materials[i].sy = vals[1];
			materials[i].sz = vals[2];
			materials[i].ex = vals[3];
			materials[i].ey = vals[4];
			materials[i].ez = vals[5];
			materials[i].ro = 0.0;
			materials[i].cv = 0.0;
		} else {
			// legacy format: sx sy sz ro [cv]
			materials[i].sx = vals[0];
			materials[i].sy = vals[1];
			materials[i].sz = vals[2];
			materials[i].ro = vals[3];
			materials[i].cv = (vals.size() == 5) ? vals[4] : 0.0;
			// set default epsilon to 1.0
			materials[i].ex = 1.0;
			materials[i].ey = 1.0;
			materials[i].ez = 1.0;
		}
	}

	int nqt = 0;
	if (idx < (int)lines.size()) {
		std::vector<int> vals;
		if (parse_ints(lines[idx], vals) && !vals.empty()) {
			nqt = vals[0];
			++idx;
		}
	}

	current_sources.assign(nodes.size(), 0.0);
	for (int i = 0; i < nqt && idx < (int)lines.size(); ++i, ++idx) {
		std::vector<double> vals;
		if (!parse_doubles(lines[idx], vals) || vals.size() < 4) break;
		int elem_id = (int)vals[0];
		double ix = vals[1], iy = vals[2], iz = vals[3];
		int e = elem_id - 1;
		if (0 <= e && e < (int)elements.size()) {
			double per_node = (ix + iy + iz) / 4.0;
			for (int k = 0; k < 4; ++k) {
				int n = elements[e].nodes[k];
				if (0 <= n && n < (int)current_sources.size()) current_sources[n] += per_node;
			}
		}
	}

	printf("Read material & BC: %d materials, %d boundary conditions, %d current sources\n", ntmate, (int)boundaries.size(), nqt);
}

void ElectrostaticAnalyzer::read_initial_potential(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	fixed_node.assign(nodes.size(), false);
	fixed_value.assign(nodes.size(), 0.0);

	int ninit = 0;
	if (fscanf(fp, "%d", &ninit) != 1) {
		fclose(fp);
		printf("Invalid header in %s\n", filename);
		return;
	}

	for (int i = 0; i < ninit; ++i) {
		int mat_id_raw = 0;
		double value = 0.0;
		if (fscanf(fp, "%d %lf", &mat_id_raw, &value) != 2) break;
		int mat_id = mat_id_raw - 1;
		fixed_materials.push_back(FixedMaterialEntry{mat_id, value});
		for (size_t e = 0; e < elements.size(); ++e) {
			if (elements[e].material == mat_id) {
				for (int k = 0; k < 4; ++k) {
					add_fixed_node(elements[e].nodes[k], value, true);
				}
			}
		}
	}

	fclose(fp);
	printf("Read initial potential: %d entries\n", ninit);
}

void ElectrostaticAnalyzer::read_analysis_config(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	std::vector<std::string> lines;
	char buf[1024];
	while (fgets(buf, sizeof(buf), fp)) {
		std::string line(buf);
		bool all_space = true;
		for (char c : line) {
			if (!isspace((unsigned char)c)) {
				all_space = false;
				break;
			}
		}
		if (!all_space) lines.push_back(line);
	}
	fclose(fp);

	auto parse_ints = [](const std::string &line, std::vector<int> &out) -> bool {
		std::istringstream iss(line);
		int v;
		out.clear();
		while (iss >> v) out.push_back(v);
		return !out.empty();
	};

	auto parse_doubles = [](const std::string &line, std::vector<double> &out) -> bool {
		std::istringstream iss(line);
		double v;
		out.clear();
		while (iss >> v) out.push_back(v);
		return !out.empty();
	};

	std::vector<int> flags;
	std::vector<int> params;
	std::vector<double> time_vals;
	std::vector<double> freq_vals;

	if (lines.size() > 0) parse_ints(lines[0], flags);
	if (lines.size() > 1) parse_ints(lines[1], params);
	if (lines.size() > 2) parse_doubles(lines[2], time_vals);
	if (lines.size() > 3) parse_doubles(lines[3], freq_vals);

	if (time_vals.size() >= 1) dt = time_vals[0];
	double time0 = (time_vals.size() >= 2) ? time_vals[1] : 0.0;
	double freq = (freq_vals.size() >= 1) ? freq_vals[0] : 0.0;

	nstep = params[0] > 0 ? params[0] : 200;
	printf("Read config: dt=%.3e, nstep=%d\n", dt, nstep);
}

void ElectrostaticAnalyzer::apply_boundary_conditions()
{
	// sin.dat/thermo.dat の Dirichlet 条件を追加
	for (size_t i = 0; i < boundaries.size(); ++i) {
		const BoundaryCondition &bc = boundaries[i];
		if (bc.type == 1 || bc.type == 4) {
			// 端点2つを固定
			add_fixed_node(bc.node1, bc.value, false);
			add_fixed_node(bc.node2, bc.value, false);
		} else {
			// その他は node1 を固定
			add_fixed_node(bc.node1, bc.value, false);
		}
	}

	for (int i = 0; i < (int)nodes.size(); ++i) {
		if (fixed_node[i]) potentials[i] = fixed_value[i];
	}
}

void ElectrostaticAnalyzer::assemble_global_matrix()
{
	auto t_start = std::chrono::high_resolution_clock::now();
	const int n = (int)nodes.size();
	global_K_rows.assign(n, {});
	global_E_rows.assign(n, {});
	global_F.assign(n, 0.0);

	for (size_t e = 0; e < elements.size(); ++e) {
		const Element &el = elements[e];
		int i0 = el.nodes[0], i1 = el.nodes[1], i2 = el.nodes[2], i3 = el.nodes[3];
		if (i0 < 0 || i1 < 0 || i2 < 0 || i3 < 0 || i0 >= n || i1 >= n || i2 >= n || i3 >= n) continue;

		const Node &n0 = nodes[i0];
		const Node &n1 = nodes[i1];
		const Node &n2 = nodes[i2];
		const Node &n3 = nodes[i3];

		double J[3][3] = {
			{n1.x - n0.x, n2.x - n0.x, n3.x - n0.x},
			{n1.y - n0.y, n2.y - n0.y, n3.y - n0.y},
			{n1.z - n0.z, n2.z - n0.z, n3.z - n0.z}
		};

		double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
					- J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
					+ J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);

		if (std::fabs(detJ) < 1e-20) continue;

		double invJ[3][3];
		invJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / detJ;
		invJ[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / detJ;
		invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ;
		invJ[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / detJ;
		invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ;
		invJ[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / detJ;
		invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ;
		invJ[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / detJ;
		invJ[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / detJ;

		// 基準四面体の dN/dxi, dN/deta, dN/dzeta
		const double dN_ref[4][3] = {
			{-1.0, -1.0, -1.0},
			{ 1.0,  0.0,  0.0},
			{ 0.0,  1.0,  0.0},
			{ 0.0,  0.0,  1.0}
		};

		// 物理座標での勾配
		double gradN[4][3];
		for (int a = 0; a < 4; ++a) {
			for (int r = 0; r < 3; ++r) {
				gradN[a][r] = invJ[r][0] * dN_ref[a][0] + invJ[r][1] * dN_ref[a][1] + invJ[r][2] * dN_ref[a][2];
			}
		}

		int mat = el.material;
		if (mat < 0 || mat >= (int)materials.size()) mat = 0;
		const Material &M = materials[mat];

		double volume = std::fabs(detJ) / 6.0;
		int gid[4] = {i0, i1, i2, i3};

		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b) {
				// anisotropic conductivity contribution
				double kxx = M.sx * gradN[a][0] * gradN[b][0];
				double kyy = M.sy * gradN[a][1] * gradN[b][1];
				double kzz = M.sz * gradN[a][2] * gradN[b][2];
				double ke_sigma = (kxx + kyy + kzz) * volume;
				global_K_rows[gid[a]].push_back(std::make_pair(gid[b], ke_sigma));

				// anisotropic permittivity (for time term matrix E)
				double e_xx = M.ex * gradN[a][0] * gradN[b][0];
				double e_yy = M.ey * gradN[a][1] * gradN[b][1];
				double e_zz = M.ez * gradN[a][2] * gradN[b][2];
				double ke_eps = (e_xx + e_yy + e_zz) * volume;
				global_E_rows[gid[a]].push_back(std::make_pair(gid[b], ke_eps));
			}
			global_F[gid[a]] += (gid[a] < (int)current_sources.size() ? current_sources[gid[a]] * volume / 4.0 : 0.0);
		}
	}

	build_sparse_matrix();

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	printf("Assembled sparse global matrix: %dx%d, nnz=%d (%.6f s)\n", n, n, (int)global_K_values.size(), elapsed);
}

void ElectrostaticAnalyzer::solve_linear_system()
{
	auto t_start = std::chrono::high_resolution_clock::now();
	const int n = (int)nodes.size();
	const bool use_ssor = (n >= 50000);//falseにすれば常にICCG、trueにすれば大規模問題でSSOR-CGを使用する
	const int max_iter = use_ssor ? std::max(1000, nstep * 10) : std::max(50, nstep);
	const double tol = use_ssor ? 1e-6 : 1e-10;
	use_ssor_preconditioner = use_ssor;

	std::vector<double> rhs = global_F;
	apply_dirichlet_constraints(rhs);

	std::vector<double> x = potentials;
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) x[i] = fixed_value[i];
	}

	if (use_ssor_preconditioner) {
		preconditioner_diag.assign(n, 0.0);
		for (int i = 0; i < n; ++i) {
			double diag = 0.0;
			for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
				if (global_K_col_idx[p] == i) {
					diag = global_K_values[p];
					break;
				}
			}
			if (!std::isfinite(diag) || std::fabs(diag) < 1e-30) diag = 1.0;
			preconditioner_diag[i] = diag;
		}
		ic_built = false;
		printf("Using SSOR preconditioner for n=%d\n", n);
	} else {
		// Build IC(0) preconditioner
		build_incomplete_cholesky();
		printf("Using IC(0) preconditioner for n=%d\n", n);
	}

	std::vector<double> r(n, 0.0), z(n, 0.0), p(n, 0.0), Ap(n, 0.0);
	matvec(x, Ap);
	for (int i = 0; i < n; ++i) {
		r[i] = rhs[i] - Ap[i];
		if (fixed_node[i]) r[i] = 0.0;
	}
	// apply preconditioner z = M^{-1} r
	apply_preconditioner(r, z);
	p = z;
	auto dot = [](const std::vector<double> &a, const std::vector<double> &b) -> double {
		double s = 0.0;
		for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
		return s;
	};

	double rnorm0 = std::sqrt(dot(r, r));
	if (rnorm0 < 1e-30) rnorm0 = 1e-30;

	double rzold = dot(r, z);
	if (std::sqrt(dot(r, r)) / rnorm0 < tol) {
		potentials = x;
		last_solver_iterations = 0;
		last_solver_method = use_ssor_preconditioner ? "SSOR-CG" : "ICCG";
		last_solver_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();
		printf("Converged at iteration 0\n");
		printf("Solved linear system with %s\n", last_solver_method.c_str());
		return;
	}

	int iter_done = 0;
	int iter_count = 0;
	for (int iter = 0; iter < max_iter; ++iter) {
		double rnorm = std::sqrt(dot(r, r));
		double rel_rnorm = rnorm / rnorm0;
		printf("ICCG iteration %d: residual=%.3e rel=%.3e\n", iter + 1, rnorm, rel_rnorm);
		matvec(p, Ap);
		double denom = dot(p, Ap);
		if (std::fabs(denom) < 1e-30) break;

		double alpha = rzold / denom;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			if (fixed_node[i]) {
				x[i] = fixed_value[i];
				r[i] = 0.0;
			}
		}

		// check convergence
		double rnorm2 = std::sqrt(dot(r, r));
		iter_done = iter + 1;
		if (rnorm2 / rnorm0 < tol) {
			printf("Converged at iteration %d\n", iter_done);
			iter_count = iter_done;
			break;
		}

		apply_preconditioner(r, z);
		double rznew = dot(r, z);
		double beta = rznew / rzold;
		for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
		rzold = rznew;
	}
	if (iter_count == 0) iter_count = iter_done;

	potentials = x;
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) potentials[i] = fixed_value[i];
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	last_solver_iterations = iter_count;
	last_solver_time = elapsed;
	last_solver_method = use_ssor_preconditioner ? "SSOR-CG" : "ICCG";
	printf("Solved linear system with %s (%.6f s)\n", last_solver_method.c_str(), elapsed);
}

void ElectrostaticAnalyzer::solve()
{
	const double conv_tol = 1e-8;
	printf("Starting electrostatic analysis...\n");
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;

	potentials.assign(nodes.size(), 0.0);
	if (fixed_node.size() != nodes.size()) {
		fixed_node.assign(nodes.size(), false);
		fixed_value.assign(nodes.size(), 0.0);
	}

	if (dt <= 0.0) {
		printf("Steady-state solve: dt=%.3e (time term disabled)\n", dt);
		apply_boundary_conditions();
		if (std::none_of(fixed_node.begin(), fixed_node.end(), [](bool v) { return v; })) {
			// 境界固定が無い場合は基準点を1つだけ 0V に固定して特異性を避ける
			if (!fixed_node.empty()) {
				fixed_node[0] = true;
				fixed_value[0] = 0.0;
				potentials[0] = 0.0;
				printf("No fixed nodes found; node 1 is used as a reference potential (0.0V)\n");
			}
		}
		assemble_global_matrix();
		solve_linear_system();
		final_step = 1;
		final_time = 0.0;
		final_max_diff = 0.0;
		step_solver_iterations.clear();
		step_solver_times.clear();
		step_solver_iterations.push_back(last_solver_iterations);
		step_solver_times.push_back(last_solver_time);
		write_potential_distribution((step_output_prefix + "001").c_str());
		printf("Analysis completed.\n");
		write_analysis_report(report_output_path.c_str());
		return;
	}

	apply_boundary_conditions();
	if (std::none_of(fixed_node.begin(), fixed_node.end(), [](bool v) { return v; })) {
		// 境界固定が無い場合は基準点を1つだけ 0V に固定して特異性を避ける
		if (!fixed_node.empty()) {
			fixed_node[0] = true;
			fixed_value[0] = 0.0;
			potentials[0] = 0.0;
			printf("No fixed nodes found; node 1 is used as a reference potential (0.0V)\n");
		}
	}
	assemble_global_matrix();
	const int n = (int)nodes.size();
	std::vector<int> perm(n, 0);
	std::vector<int> inv_perm(n, 0);
	if (n > 0) {
		// Reverse Cuthill-McKee ordering on the assembled adjacency.
		std::vector<int> degree(n, 0);
		for (int i = 0; i < n; ++i) {
			int deg = 0;
			for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
				if (global_K_col_idx[p] != i) ++deg;
			}
			degree[i] = deg;
		}

		std::vector<char> visited(n, 0);
		std::vector<int> rcm_order;
		rcm_order.reserve(n);
		std::vector<int> frontier;
		frontier.reserve(n);
		std::vector<int> neighbors;
		neighbors.reserve(64);

		while ((int)rcm_order.size() < n) {
			int start = -1;
			int best_degree = 0;
			for (int i = 0; i < n; ++i) {
				if (visited[i]) continue;
				if (start < 0 || degree[i] < best_degree) {
					start = i;
					best_degree = degree[i];
				}
			}
			if (start < 0) break;

			frontier.clear();
			frontier.push_back(start);
			visited[start] = 1;

			for (size_t head = 0; head < frontier.size(); ++head) {
				int node = frontier[head];
				rcm_order.push_back(node);

				neighbors.clear();
				for (int p = global_K_row_ptr[node]; p < global_K_row_ptr[node + 1]; ++p) {
					int nb = global_K_col_idx[p];
					if (nb == node || visited[nb]) continue;
					neighbors.push_back(nb);
					visited[nb] = 1;
				}
				std::sort(neighbors.begin(), neighbors.end(), [&](int lhs, int rhs) {
					if (degree[lhs] != degree[rhs]) return degree[lhs] < degree[rhs];
					return lhs < rhs;
				});
				for (int nb : neighbors) {
					frontier.push_back(nb);
				}
			}
		}

		std::reverse(rcm_order.begin(), rcm_order.end());
		for (int i = 0; i < n; ++i) {
			perm[i] = rcm_order[i];
			inv_perm[perm[i]] = i;
		}
	}

	auto permute_vector = [&](const std::vector<double> &src, std::vector<double> &dst) {
		dst.assign(n, 0.0);
		for (int i = 0; i < n; ++i) {
			dst[i] = src[perm[i]];
		}
	};

	auto permute_bool_vector = [&](const std::vector<bool> &src, std::vector<bool> &dst) {
		dst.assign(n, false);
		for (int i = 0; i < n; ++i) {
			dst[i] = src[perm[i]];
		}
	};

	auto permute_csr = [&](const std::vector<int> &src_row_ptr, const std::vector<int> &src_col_idx,
		const std::vector<double> &src_values,
		std::vector<int> &dst_row_ptr, std::vector<int> &dst_col_idx, std::vector<double> &dst_values) {
		dst_row_ptr.assign(n + 1, 0);
		dst_col_idx.clear();
		dst_values.clear();
		std::vector<std::pair<int, double>> row;
		for (int new_row = 0; new_row < n; ++new_row) {
			int old_row = perm[new_row];
			row.clear();
			for (int p = src_row_ptr[old_row]; p < src_row_ptr[old_row + 1]; ++p) {
				int old_col = src_col_idx[p];
				row.push_back(std::make_pair(inv_perm[old_col], src_values[p]));
			}
			std::sort(row.begin(), row.end(), [](const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) {
				return lhs.first < rhs.first;
			});
			for (size_t i = 0; i < row.size(); ) {
				int col = row[i].first;
				double sum = row[i].second;
				size_t j = i + 1;
				while (j < row.size() && row[j].first == col) {
					sum += row[j].second;
					++j;
				}
				if (std::fabs(sum) > 1e-30) {
					dst_col_idx.push_back(col);
					dst_values.push_back(sum);
				}
				i = j;
			}
			dst_row_ptr[new_row + 1] = (int)dst_col_idx.size();
		}
	};

	const std::vector<int> original_K_row_ptr = global_K_row_ptr;
	const std::vector<int> original_K_col_idx = global_K_col_idx;
	const std::vector<double> original_K_values = global_K_values;
	const std::vector<int> original_E_row_ptr = global_E_row_ptr;
	const std::vector<int> original_E_col_idx = global_E_col_idx;
	const std::vector<double> original_E_values = global_E_values;

	permute_csr(original_K_row_ptr, original_K_col_idx, original_K_values, global_K_row_ptr, global_K_col_idx, global_K_values);
	permute_csr(original_E_row_ptr, original_E_col_idx, original_E_values, global_E_row_ptr, global_E_col_idx, global_E_values);

	std::vector<double> permuted_F;
	permute_vector(global_F, permuted_F);
	global_F.swap(permuted_F);

	std::vector<double> permuted_potentials;
	permute_vector(potentials, permuted_potentials);
	potentials.swap(permuted_potentials);

	std::vector<bool> permuted_fixed_node;
	permute_bool_vector(fixed_node, permuted_fixed_node);
	fixed_node.swap(permuted_fixed_node);

	std::vector<double> permuted_fixed_value;
	permute_vector(fixed_value, permuted_fixed_value);
	fixed_value.swap(permuted_fixed_value);

	std::vector<double> base_F = global_F;
	printf("Transient solve: dt=%.3e, nstep=%d\n", dt, nstep);
	step_solver_iterations.clear();
	step_solver_times.clear();

	std::vector<double> potentials_prev = potentials;
	for (int step = 0; step < nstep; ++step) {
		std::vector<double> potentials_old = potentials_prev;
		std::vector<int> step_row_ptr;
		std::vector<int> step_col_idx;
		std::vector<double> step_values;
		build_time_step_matrix(1.0 / dt, step_row_ptr, step_col_idx, step_values);

		global_K_row_ptr.swap(step_row_ptr);
		global_K_col_idx.swap(step_col_idx);
		global_K_values.swap(step_values);

		global_F = base_F;
		std::vector<double> eps_phi;
		matvec_csr(global_E_row_ptr, global_E_col_idx, global_E_values, potentials_prev, eps_phi);
		for (size_t i = 0; i < global_F.size(); ++i) {
			global_F[i] += eps_phi[i] / dt;
		}

		potentials = potentials_prev;
		solve_linear_system();
		potentials_prev = potentials;
		step_solver_iterations.push_back(last_solver_iterations);
		step_solver_times.push_back(last_solver_time);

		double max_diff = 0.0;
		for (size_t i = 0; i < potentials.size(); ++i) {
			double diff = std::fabs(potentials[i] - potentials_old[i]);
			if (diff > max_diff) max_diff = diff;
		}

		double vmin = 0.0;
		double vmax = 0.0;
		if (!potentials.empty()) {
			vmin = *std::min_element(potentials.begin(), potentials.end());
			vmax = *std::max_element(potentials.begin(), potentials.end());
		}
		printf("Step %d/%d: t=%.6e s, Vmin=%.6e, Vmax=%.6e, max_diff=%.6e\n",
			step + 1, nstep, (step + 1) * dt, vmin, vmax, max_diff);

		char outname[64];
		std::snprintf(outname, sizeof(outname), "%s%03d", step_output_prefix.c_str(), step + 1);
		write_potential_distribution(outname);

		final_step = step + 1;
		final_time = (step + 1) * dt;
		final_max_diff = max_diff;

		if (max_diff < conv_tol) {
			printf("Converged at step %d: max_diff=%.6e < tol=%.6e\n", step + 1, max_diff, conv_tol);
			converged_early = true;
			break;
		}
	}

	// Restore original node ordering for output/reporting.
	std::vector<double> original_potentials(n, 0.0);
	for (int i = 0; i < n; ++i) {
		original_potentials[perm[i]] = potentials[i];
	}
	potentials.swap(original_potentials);

	printf("Analysis completed.\n");
	write_analysis_report(report_output_path.c_str());
}

void ElectrostaticAnalyzer::write_analysis_report(const char *filename) const
{
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	fprintf(fp, "Electrostatic FEM Analysis Summary\n");
	fprintf(fp, "=================================\n\n");
	fprintf(fp, "Mesh\n");
	fprintf(fp, "----\n");
	fprintf(fp, "Nodes: %d\n", (int)nodes.size());
	fprintf(fp, "Elements: %d\n\n", (int)elements.size());

	fprintf(fp, "Analysis Conditions\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "dt: %.6e\n", dt);
	fprintf(fp, "nstep(max): %d\n", nstep);
	fprintf(fp, "convergence_tolerance: %.6e\n", 1e-8);
	fprintf(fp, "final_step: %d\n", final_step);
	fprintf(fp, "final_time: %.6e\n", final_time);
	fprintf(fp, "final_max_diff: %.6e\n", final_max_diff);
	fprintf(fp, "converged_early: %s\n\n", converged_early ? "yes" : "no");
	fprintf(fp, "solver_method: %s\n", last_solver_method.c_str());
	fprintf(fp, "last_solver_iterations: %d\n", last_solver_iterations);
	fprintf(fp, "last_solver_time_sec: %.6e\n\n", last_solver_time);

	fprintf(fp, "Solver Iterations Per Step\n");
	fprintf(fp, "-------------------------\n");
	if (step_solver_iterations.empty()) {
		fprintf(fp, "(no iteration history)\n\n");
	} else {
		fprintf(fp, "step,iterations,solve_time_sec\n");
		for (size_t i = 0; i < step_solver_iterations.size(); ++i) {
			double solve_time = (i < step_solver_times.size()) ? step_solver_times[i] : 0.0;
			fprintf(fp, "%zu,%d,%.6e\n", i + 1, step_solver_iterations[i], solve_time);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "Boundary Conditions\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "count: %d\n", (int)boundaries.size());
	for (size_t i = 0; i < boundaries.size(); ++i) {
		const BoundaryCondition &bc = boundaries[i];
		fprintf(fp, "%zu: node1=%d node2=%d type=%d value=%.6e\n",
			i + 1, bc.node1 + 1, bc.node2 + 1, bc.type, bc.value);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Fixed Materials (tmate.dat)\n");
	fprintf(fp, "---------------------------\n");
	fprintf(fp, "count: %d\n", (int)fixed_materials.size());
	for (size_t i = 0; i < fixed_materials.size(); ++i) {
		const FixedMaterialEntry &fm = fixed_materials[i];
		fprintf(fp, "%zu: material_id=%d value=%.6e\n", i + 1, fm.material_id + 1, fm.value);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Material Properties\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "count: %d\n", (int)materials.size());
	for (size_t i = 0; i < materials.size(); ++i) {
		const Material &m = materials[i];
		fprintf(fp, "%zu: sigma=(%.6e, %.6e, %.6e) epsilon=(%.6e, %.6e, %.6e)\n",
			i + 1, m.sx, m.sy, m.sz, m.ex, m.ey, m.ez);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Material Block Ranges\n");
	fprintf(fp, "---------------------\n");
	fprintf(fp, "count: %d\n", (int)material_ranges.size());
	for (size_t i = 0; i < material_ranges.size(); ++i) {
		const MaterialBlockRange &r = material_ranges[i];
		fprintf(fp, "%zu: start=%d end=%d material_id=%d\n",
			i + 1, r.start_element, r.end_element, r.material_id + 1);
	}
	fclose(fp);
	printf("Wrote %s\n", filename);
}

void ElectrostaticAnalyzer::write_potential_distribution(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	auto format_e12_4 = [](double value) -> std::string {
		char buf[64];
		char out[64];

		if (!std::isfinite(value)) {
			value = 0.0;
		}

		// Keep the formatted field within 12 characters.
		// Values that would need a 3-digit exponent are treated as signed zero.
		if (std::fabs(value) < 1.0e-99) {
			std::snprintf(buf, sizeof(buf), "%s0.0000E+00", std::signbit(value) ? "-" : " ");
		} else {
			std::snprintf(buf, sizeof(buf), "%.4E", value);
			char *e_pos = std::strchr(buf, 'E');
			if (e_pos != nullptr) {
				int exponent = std::atoi(e_pos + 1);
				if (std::abs(exponent) > 99) {
					std::snprintf(buf, sizeof(buf), "%s0.0000E+00", std::signbit(value) ? "-" : " ");
				}
			}
		}

		std::snprintf(out, sizeof(out), "%12s", buf);
		return std::string(out);
	};

	int col = 0;
	for (size_t i = 0; i < potentials.size(); ++i) {
		fprintf(fp, "%s", format_e12_4(potentials[i]).c_str());
		col++;
		if (col == 6 || i + 1 == potentials.size()) {
			fprintf(fp, "\n");
			col = 0;
		}
	}

	fclose(fp);
	printf("Wrote %s: %d nodal potentials\n", filename, (int)potentials.size());
}
