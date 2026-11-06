#include "electrostatic_analyzer.h"

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <vector>
#include <map>
#include <set>
#include <array>
#include <limits>

// ============================================================================
// 初期化・共通設定
// ============================================================================

ElectrostaticAnalyzer::ElectrostaticAnalyzer()
	: mesh_node_count(0), mesh_element_count(0), dt(0.0), nstep(200), unit(1.0)
{
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;
	last_solver_iterations = 0;
	last_solver_time = 0.0;
	last_solver_relative_residual = 0.0;
	last_solver_converged = false;
	last_solver_method = "ICCG";
	step_output_prefix = "tempa.dat";
	report_output_path = "analysis_summary.txt";
}

void ElectrostaticAnalyzer::set_output_paths(const std::string &step_prefix, const std::string &report_path)
{
	step_output_prefix = step_prefix.empty() ? std::string("tempa.dat") : step_prefix;
	report_output_path = report_path.empty() ? std::string("analysis_summary.txt") : report_path;
}

void ElectrostaticAnalyzer::read_voltage_difference_config(const char *filename)
{
	std::ifstream input(filename);
	if (!input) {
		printf("Cannot open %s\n", filename);
		return;
	}

	voltage_difference_monitors.clear();
	std::string line;
	while (std::getline(input, line)) {
		size_t first = line.find_first_not_of(" \t\r\n");
		if (first == std::string::npos || line[first] == '#') continue;

		std::istringstream row(line);
		int node_a = 0;
		int node_b = 0;
		std::string label;
		if (!(row >> node_a >> node_b)) continue;
		std::getline(row, label);
		size_t label_start = label.find_first_not_of(" \t");
		label = label_start == std::string::npos
			? "node_" + std::to_string(node_b) + "_minus_node_" + std::to_string(node_a)
			: label.substr(label_start);

		if (node_a < 1 || node_a > mesh_node_count || node_b < 1 || node_b > mesh_node_count) {
			printf("Ignoring voltage difference monitor with invalid nodes: %d %d\n", node_a, node_b);
			continue;
		}
		voltage_difference_monitors.push_back(VoltageDifferenceMonitor{node_a - 1, node_b - 1, label, {}, {}});
	}
	printf("Read %d voltage difference monitor(s) from %s\n",
		(int)voltage_difference_monitors.size(), filename);
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

	// global_Kとglobal_Eは、各行の列番号順に整列・集約済みである。
	// 行コピーと再ソートを避けるため、2つのCSR行列を直接マージする。
	// scale_eps=1/dtの場合、A = K_sigma + K_epsilon/dt を作る。
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
	// Assembly contributions are no longer needed after CSR conversion.
	decltype(global_K_rows)().swap(global_K_rows);
	decltype(global_E_rows)().swap(global_E_rows);

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	printf("build_sparse_matrix: %.6f s\n", elapsed);
}

void ElectrostaticAnalyzer::apply_dirichlet_constraints(std::vector<double> &rhs)
{
	const int n = (int)nodes.size();
	// 固定節点jの既知電位を、未知節点iの右辺へ移す：
	//   rhs[i] <- rhs[i] - A[i,j] * fixed_value[j]
	// 同時に固定節点列の非対角成分を0にする。計算量はO(nnz)。
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

	// 固定節点の式を phi[j] = fixed_value[j] に置き換える。
	// そのため固定節点行は、対角成分を1、その他を0とする。
	// 行と列の両方を処理することで、CG法に必要な行列の対称性を保つ。
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

	// IC(0)分解前に対称対角スケーリングを行い、係数の桁差を小さくする。
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

	std::vector<int> source_positions;
	std::vector<int> diagonal_positions(n, -1);
	// 下三角部分（対角を含む、列番号<=行番号）の非ゼロパターンを作る。
	for (int i = 0; i < n; ++i) {
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			int c = global_K_col_idx[p];
			if (c <= i) {
				source_positions.push_back(p);
				if (c == i) diagonal_positions[i] = (int)L_col_idx.size();
				L_col_idx.push_back(c);
				// 数値は後のIC(0)計算で設定するため、ここでは仮に0を入れる。
				L_values.push_back(0.0);
			}
		}
		L_row_ptr[i + 1] = (int)L_col_idx.size();
	}

	// A ≒ L L^T となるLを、元行列と同じ非ゼロパターンの範囲で計算する。
	for (int i = 0; i < n; ++i) {
		// 非対角成分L(i,j)、j<iを計算する。
		for (int p = L_row_ptr[i]; p < L_row_ptr[i+1]; ++p) {
			int j = L_col_idx[p];
			if (j == i) continue;
			double val = ic_scale[i] * global_K_values[source_positions[p]] * ic_scale[j];

			// jより前にある共通の非ゼロ列について、L(i,k)L(j,k)を差し引く。
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

			// 行jの対角成分L(j,j)を探す。
			double Ljj = 0.0;
			if (diagonal_positions[j] >= 0) Ljj = L_values[diagonal_positions[j]];
			if (std::fabs(Ljj) < 1e-30) Ljj = 1e-30;
			L_values[p] = val / Ljj;
		}

		// 対角成分L(i,i)を計算する。
		const int diagonal = diagonal_positions[i];
		double diag = diagonal >= 0 ? ic_scale[i] * global_K_values[source_positions[diagonal]] * ic_scale[i] : 0.0;
		for (int q = L_row_ptr[i]; q < L_row_ptr[i+1]; ++q) {
			int k = L_col_idx[q];
			if (k == i) continue;
			diag -= L_values[q] * L_values[q];
		}
		if (!std::isfinite(diag) || diag <= 0.0) diag = 1e-30;
		// 計算した対角成分を格納する。
		if (diagonal >= 0) L_values[diagonal] = std::sqrt(diag);
	}

	// L^Tによる後退代入を効率化するため、LのCSC（列圧縮）形式も作る。
	L_csc_col_ptr.assign(n + 1, 0);
	for (size_t idx = 0; idx < L_col_idx.size(); ++idx) {
		int col = L_col_idx[idx];
		++L_csc_col_ptr[col + 1];
	}
	for (int i = 0; i < n; ++i) L_csc_col_ptr[i+1] += L_csc_col_ptr[i];
	L_csc_row_idx.assign(L_col_idx.size(), 0);
	L_csc_values.assign(L_values.size(), 0.0);
	// CSR形式のLをCSC形式へ詰め替える。
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

	// 前進スイープ：(D + L)y = r を解く。
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

	// 中間処理として対角行列Dを掛ける。
	for (int i = 0; i < n; ++i) {
		y[i] *= diag_values[i];
	}

	// 後退スイープ：(D + U)z = Dy を解く。
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
		// IC(0)を構築できていない場合は、対角前処理へ切り替える。
		for (int i = 0; i < n; ++i) {
			double Aii = 0.0;
			for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i+1]; ++p) if (global_K_col_idx[p] == i) { Aii = global_K_values[p]; break; }
			if (std::fabs(Aii) < 1e-30) z[i] = r[i]; else z[i] = r[i] / Aii;
		}
		return;
	}

	// Reuse scratch storage; forward substitution overwrites each entry before use.
	preconditioner_work.resize(n);
	auto &y = preconditioner_work;

	// 前進代入で Ly = r を解く。
	for (int i = 0; i < n; ++i) {
		double s = r[i] * ic_scale[i];
		double Lii = 1.0;
		for (int p = L_row_ptr[i]; p < L_row_ptr[i+1]; ++p) {
			int c = L_col_idx[p];
			if (c == i) { Lii = L_values[p]; break; }
			s -= L_values[p] * y[c];
		}
		if (std::fabs(Lii) < 1e-30) Lii = 1e-30;
		y[i] = s / Lii;
	}

	// CSC形式を使った後退代入で L^T z = y を解く。
	for (int col = n - 1; col >= 0; --col) {
		double s = y[col];
		double Lii = 1.0;
		// CSCの列colから、row==colとなる対角成分を探す。
		for (int p = L_csc_col_ptr[col]; p < L_csc_col_ptr[col+1]; ++p) {
			int row = L_csc_row_idx[p];
			if (row == col) { Lii = L_csc_values[p]; break; }
		}
		for (int p = L_csc_col_ptr[col]; p < L_csc_col_ptr[col+1]; ++p) {
			int row = L_csc_row_idx[p];
			if (row == col) continue;
			// L(row,col) * z[row]を現在の式から差し引く。
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
	// 外部公開用のgetterが返す節点数・要素数を保存する。
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

	// 材料値を初期化する。既定値は導電率・誘電率ともに1とする。
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
			// 旧形式：sx sy sz ro [cv]
			materials[i].sx = vals[0];
			materials[i].sy = vals[1];
			materials[i].sz = vals[2];
			materials[i].ro = vals[3];
			materials[i].cv = (vals.size() == 5) ? vals[4] : 0.0;
			// 旧形式には誘電率がないため、既定値1.0を設定する。
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

void ElectrostaticAnalyzer::read_pulse_config(const char *filename)
{
	FILE *fp = fopen(filename, "r");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	pulse_voltages.clear();
	char buf[1024];
	int line_no = 0;
	while (fgets(buf, sizeof(buf), fp)) {
		++line_no;
		std::string line(buf);
		size_t comment = line.find('#');
		if (comment != std::string::npos) line = line.substr(0, comment);
		bool all_space = true;
		for (char c : line) {
			if (!isspace((unsigned char)c)) {
				all_space = false;
				break;
			}
		}
		if (all_space) continue;

		int mat_id_raw = 0;
		PulseVoltageEntry pulse;
		pulse.phase = 0.0;
		int nread = std::sscanf(line.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf",
			&mat_id_raw, &pulse.v_low, &pulse.v_high, &pulse.t_rise,
			&pulse.t_high, &pulse.t_fall, &pulse.t_low, &pulse.phase);
		if (nread != 7 && nread != 8) {
			printf("Invalid pulse record at line %d in %s\n", line_no, filename);
			continue;
		}
		pulse.material_id = mat_id_raw - 1;
		if (pulse.material_id < 0) {
			printf("Invalid pulse material id at line %d in %s\n", line_no, filename);
			continue;
		}
		if (pulse.t_rise < 0.0 || pulse.t_high < 0.0 || pulse.t_fall < 0.0 || pulse.t_low < 0.0) {
			printf("Invalid negative pulse duration at line %d in %s\n", line_no, filename);
			continue;
		}

		for (size_t e = 0; e < elements.size(); ++e) {
			if (elements[e].material == pulse.material_id) {
				for (int k = 0; k < 4; ++k) {
					int node = elements[e].nodes[k];
					if (0 <= node && node < (int)nodes.size()) {
						pulse.nodes.push_back(node);
						add_fixed_node(node, pulse.v_low, true);
					}
				}
			}
		}
		std::sort(pulse.nodes.begin(), pulse.nodes.end());
		pulse.nodes.erase(std::unique(pulse.nodes.begin(), pulse.nodes.end()), pulse.nodes.end());
		pulse_voltages.push_back(pulse);
	}

	fclose(fp);
	printf("Read pulse config: %d entries\n", (int)pulse_voltages.size());
}

double ElectrostaticAnalyzer::evaluate_pulse_voltage(const PulseVoltageEntry &pulse, double time) const
{
	const double period = pulse.t_rise + pulse.t_high + pulse.t_fall + pulse.t_low;
	if (period <= 0.0) return pulse.v_low;

	double t = std::fmod(time - pulse.phase, period);
	if (t < 0.0) t += period;

	if (pulse.t_rise > 0.0 && t < pulse.t_rise) {
		double a = t / pulse.t_rise;
		return pulse.v_low + (pulse.v_high - pulse.v_low) * a;
	}
	t -= pulse.t_rise;

	if (t < pulse.t_high) return pulse.v_high;
	t -= pulse.t_high;

	if (pulse.t_fall > 0.0 && t < pulse.t_fall) {
		double a = t / pulse.t_fall;
		return pulse.v_high + (pulse.v_low - pulse.v_high) * a;
	}

	return pulse.v_low;
}

void ElectrostaticAnalyzer::apply_pulse_voltages(double time)
{
	for (const PulseVoltageEntry &pulse : pulse_voltages) {
		if (circuit.enabled && pulse.material_id == circuit.material_id) continue;
		double value = evaluate_pulse_voltage(pulse, time);
		for (int node : pulse.nodes) {
			if (0 <= node && node < (int)fixed_value.size()) {
				fixed_node[node] = true;
				fixed_value[node] = value;
				if (node < (int)potentials.size()) potentials[node] = value;
			}
		}
	}
}

namespace {
[[noreturn]] void electrical_error(const std::string &message)
{
	fprintf(stderr, "Electrical configuration error: %s\n", message.c_str());
	std::exit(EXIT_FAILURE);
}

std::map<std::string, std::string> electrical_options(const char *filename)
{
	std::ifstream file(filename);
	if (!file) electrical_error(std::string("cannot read ") + filename);
	std::map<std::string, std::string> result;
	std::string line, key, value, extra;
	while (std::getline(file, line)) {
		line = line.substr(0, line.find('#'));
		std::istringstream row(line);
		if (!(row >> key)) continue;
		if (!(row >> value) || (row >> extra) || result.count(key))
			electrical_error(std::string("invalid or duplicate option in ") + filename + ": " + key);
		result[key] = value;
	}
	return result;
}

double electrical_number(const std::string &text)
{
	char *end = nullptr;
	double value = std::strtod(text.c_str(), &end);
	if (end == text.c_str() || *end != '\0' || !std::isfinite(value))
		electrical_error("invalid finite number: " + text);
	return value;
}

void add_csr_diagonal(int row, double value, std::vector<int> &rp,
	std::vector<int> &ci, std::vector<double> &values)
{
	if (value == 0.0) return;
	int p = rp[row];
	while (p < rp[row + 1] && ci[p] < row) ++p;
	if (p < rp[row + 1] && ci[p] == row) values[p] += value;
	else {
		ci.insert(ci.begin() + p, row);
		values.insert(values.begin() + p, value);
		for (size_t r = row + 1; r < rp.size(); ++r) ++rp[r];
	}
}
}

void ElectrostaticAnalyzer::read_symmetry_config(const char *filename)
{
	auto options = electrical_options(filename);
	if (options.size() != 1 || options["axis"] != "z")
		electrical_error("symmetry.dat must contain only: axis z");
	axial_symmetry = true;
}

void ElectrostaticAnalyzer::read_circuit_config(const char *filename)
{
	auto options = electrical_options(filename);
	const std::set<std::string> allowed = {"electrode_material", "series_resistance_ohm",
		"off_resistance_ohm", "disconnect_time_s", "external_capacitance_F",
		"coordinate_scale_to_m", "initial_state"};
	for (const auto &option : options) if (!allowed.count(option.first))
		electrical_error("unknown circuit option: " + option.first);
	for (const auto &key : {"electrode_material", "series_resistance_ohm", "disconnect_time_s",
		"coordinate_scale_to_m", "initial_state"}) if (!options.count(key))
		electrical_error(std::string("missing circuit option: ") + key);
	double material = electrical_number(options["electrode_material"]);
	if (material < 1 || material > materials.size() || material != std::floor(material))
		electrical_error("invalid circuit electrode material");
	circuit.material_id = (int)material - 1;
	auto resistance = [&](const std::string &value) {
		if (value == "open") return std::numeric_limits<double>::infinity();
		if (value == "ideal") return 0.0;
		double r = electrical_number(value);
		if (r < 0) electrical_error("resistance must be nonnegative or open");
		return r;
	};
	circuit.series_resistance = resistance(options["series_resistance_ohm"]);
	circuit.disconnect_time = options["disconnect_time_s"] == "never"
		? std::numeric_limits<double>::infinity() : electrical_number(options["disconnect_time_s"]);
	if (std::isfinite(circuit.disconnect_time) && !options.count("off_resistance_ohm"))
		electrical_error("off_resistance_ohm is required when disconnecting");
	circuit.off_resistance = options.count("off_resistance_ohm")
		? resistance(options["off_resistance_ohm"]) : std::numeric_limits<double>::infinity();
	circuit.external_capacitance = options.count("external_capacitance_F")
		? electrical_number(options["external_capacitance_F"]) : 0.0;
	circuit.coordinate_scale = electrical_number(options["coordinate_scale_to_m"]);
	if (circuit.coordinate_scale <= 0 || circuit.external_capacitance < 0)
		electrical_error("coordinate scale must be positive and capacitance nonnegative");
	if (options["initial_state"] != "steady" && options["initial_state"] != "charge_free")
		electrical_error("initial_state must be steady or charge_free");
	circuit.steady_initial_state = options["initial_state"] == "steady";
	circuit.nodes.clear();
	for (const auto &el : elements) if (el.material == circuit.material_id)
		for (int node : el.nodes) {
			if (node < 0 || node >= (int)nodes.size()) electrical_error("invalid circuit electrode connectivity");
			circuit.nodes.push_back(node);
		}
	std::sort(circuit.nodes.begin(), circuit.nodes.end());
	circuit.nodes.erase(std::unique(circuit.nodes.begin(), circuit.nodes.end()), circuit.nodes.end());
	int count = 0;
	for (const auto &pulse : pulse_voltages) if (pulse.material_id == circuit.material_id) ++count;
	if (circuit.nodes.empty() || count != 1)
		electrical_error("circuit electrode needs an existing material and exactly one pulse.dat entry");
	for (const auto &bc : boundaries) if (bc.type == 1 || bc.type == 4)
		if (std::binary_search(circuit.nodes.begin(), circuit.nodes.end(), bc.node1) ||
			std::binary_search(circuit.nodes.begin(), circuit.nodes.end(), bc.node2))
			electrical_error("explicit fixed boundary conflicts with the circuit electrode");
	for (const auto &fixed : fixed_materials) if (fixed.material_id != circuit.material_id)
		for (const auto &el : elements) if (el.material == fixed.material_id)
			for (int node : el.nodes) if (std::binary_search(circuit.nodes.begin(), circuit.nodes.end(), node))
				electrical_error("circuit electrode shares a node with another fixed material");
	for (double source : current_sources) if (source != 0.0)
		electrical_error("nonzero distributed current sources are not supported with circuit coupling");
	circuit.enabled = true;
}

double ElectrostaticAnalyzer::circuit_source_voltage(double time) const
{
	if (time >= circuit.disconnect_time) return 0.0;
	for (const auto &pulse : pulse_voltages)
		if (pulse.material_id == circuit.material_id) return evaluate_pulse_voltage(pulse, time);
	electrical_error("missing circuit source waveform");
}

double ElectrostaticAnalyzer::circuit_step_time(double time) const
{
	if (!circuit.enabled) return time;
	// At the switch instant the capacitor voltage is continuous: integrate the
	// preceding interval with its ON topology, then use OFF on the next interval.
	double tolerance = 64 * std::numeric_limits<double>::epsilon()
		* std::max(std::fabs(time), dt);
	if (std::isfinite(circuit.disconnect_time) && std::fabs(time - circuit.disconnect_time) <= tolerance)
		return circuit.disconnect_time - 2 * tolerance;
	return time;
}

void ElectrostaticAnalyzer::apply_circuit_boundary(double time)
{
	if (!circuit.enabled) return;
	double resistance = time < circuit.disconnect_time ? circuit.series_resistance : circuit.off_resistance;
	double value = circuit_source_voltage(time);
	auto apply = [&](int node) {
		fixed_node[node] = resistance == 0.0;
		if (fixed_node[node]) { fixed_value[node] = value; potentials[node] = value; }
	};
	if (circuit.dof >= 0) apply(circuit.dof);
	else for (int node : circuit.nodes) apply(node);
}

void ElectrostaticAnalyzer::add_circuit_terms(double time)
{
	if (!circuit.enabled) return;
	double resistance = time < circuit.disconnect_time ? circuit.series_resistance : circuit.off_resistance;
	if (resistance > 0.0 && std::isfinite(resistance)) {
		double conductance = 1.0 / resistance;
		if (!std::isfinite(conductance)) electrical_error("resistance is too small");
		add_csr_diagonal(circuit.dof, conductance, global_K_row_ptr, global_K_col_idx, global_K_values);
		global_F[circuit.dof] += conductance * circuit_source_voltage(time);
	}
}

int ElectrostaticAnalyzer::monitor_dof(int node) const
{
	return output_node_to_dof.empty() ? node : output_node_to_dof[node];
}

void ElectrostaticAnalyzer::prepare_reduced_system()
{
	last_dof_count = (int)nodes.size();
	if (!axial_symmetry && !circuit.enabled) return;
	global_K_rows.clear(); global_K_rows.shrink_to_fit();
	global_E_rows.clear(); global_E_rows.shrink_to_fit();
	const int original_count = (int)nodes.size();
	for (const auto &el : elements) {
		if (el.material < 0 || el.material >= (int)materials.size()) electrical_error("invalid material in reduced model");
		for (int node : el.nodes) if (node < 0 || node >= original_count)
			electrical_error("invalid connectivity in reduced model");
	}
	std::vector<int> parent(original_count);
	for (int i = 0; i < original_count; ++i) parent[i] = i;
	auto find = [&](int node) {
		int root = node;
		while (parent[root] != root) root = parent[root];
		while (parent[node] != node) { int next = parent[node]; parent[node] = root; node = next; }
		return root;
	};
	if (axial_symmetry) {
		std::map<std::pair<double, double>, int> columns;
		std::set<double> levels;
		std::vector<std::vector<int>> signatures(original_count);
		for (const auto &el : elements) for (int node : el.nodes) signatures[node].push_back(el.material);
		for (auto &signature : signatures) {
			std::sort(signature.begin(), signature.end());
			signature.erase(std::unique(signature.begin(), signature.end()), signature.end());
		}
		std::vector<int> counts(original_count, 0);
		for (int i = 0; i < original_count; ++i) {
			const auto &p = nodes[i];
			if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
				electrical_error("nonfinite symmetry coordinates");
			levels.insert(p.z);
			auto result = columns.emplace(std::make_pair(p.x, p.y), i);
			int representative = result.first->second;
			if (fixed_node[i] != fixed_node[representative] ||
				(fixed_node[i] && fixed_value[i] != fixed_value[representative]) ||
				signatures[i] != signatures[representative])
				electrical_error("z symmetry conflicts with material or voltage boundary conditions");
			if (current_sources.size() == nodes.size() && current_sources[i] != current_sources[representative])
				electrical_error("z symmetry conflicts with current sources");
			parent[i] = representative;
			++counts[representative];
		}
		if (levels.size() < 2) electrical_error("z symmetry needs at least two layers");
		for (const auto &column : columns) if (counts[column.second] != (int)levels.size())
			electrical_error("z symmetry requires identical XY coordinates on every layer");
		std::map<double, int> level_index;
		for (double z : levels) level_index[z] = (int)level_index.size();
		std::map<std::array<int, 3>, int> triangle_materials;
		for (const auto &el : elements) {
			std::array<int, 4> projected;
			int low = (int)levels.size(), high = -1;
			for (int j = 0; j < 4; ++j) {
				projected[j] = parent[el.nodes[j]];
				int layer = level_index[nodes[el.nodes[j]].z];
				low = std::min(low, layer); high = std::max(high, layer);
			}
			std::sort(projected.begin(), projected.end());
			if (std::unique(projected.begin(), projected.end()) - projected.begin() != 3 || high - low != 1)
				electrical_error("z symmetry requires a conforming extruded triangular mesh");
			std::array<int, 3> triangle = {projected[0], projected[1], projected[2]};
			auto result = triangle_materials.emplace(triangle, el.material);
			if (!result.second && result.first->second != el.material)
				electrical_error("material varies along the extrusion axis");
		}
	}
	if (circuit.enabled) {
		int root = find(circuit.nodes.front());
		for (int node : circuit.nodes) parent[find(node)] = root;
	}
	output_node_to_dof.resize(original_count);
	std::map<int, int> groups;
	std::vector<Node> reduced_nodes;
	std::vector<bool> reduced_fixed;
	std::vector<double> reduced_values, reduced_potentials, reduced_F;
	for (int i = 0; i < original_count; ++i) {
		int root = find(i);
		auto result = groups.emplace(root, (int)groups.size());
		int group = result.first->second;
		output_node_to_dof[i] = group;
		if (result.second) {
			reduced_nodes.push_back(nodes[i]); reduced_fixed.push_back(false);
			reduced_values.push_back(0); reduced_potentials.push_back(potentials[i]); reduced_F.push_back(0);
		}
		if (fixed_node[i]) {
			if (reduced_fixed[group] && reduced_values[group] != fixed_value[i])
				electrical_error("incompatible fixed voltages in an equipotential group");
			reduced_fixed[group] = true; reduced_values[group] = fixed_value[i];
			reduced_potentials[group] = fixed_value[i];
		}
		reduced_F[group] += global_F[i];
	}
	input_nodes.swap(nodes);
	nodes.swap(reduced_nodes);
	auto reduce_matrix = [&](std::vector<int> &rp, std::vector<int> &ci, std::vector<double> &values) {
		std::vector<std::vector<std::pair<int, double>>> rows(nodes.size());
		for (int i = 0; i < original_count; ++i) for (int p = rp[i]; p < rp[i + 1]; ++p)
			rows[output_node_to_dof[i]].push_back({output_node_to_dof[ci[p]], values[p]});
		build_sparse_from_rows(rows, rp, ci, values);
		// FEM gradient matrices annihilate constants. Reconstruct each diagonal
		// from its off-diagonal entries so that merging a highly conductive
		// electrode does not leave a cancellation error as a spurious shunt.
		for (int row = 0; row < (int)nodes.size(); ++row) {
			long double sum = 0;
			int diagonal = -1;
			for (int p = rp[row]; p < rp[row + 1]; ++p) {
				if (ci[p] == row) diagonal = p;
				else sum += values[p];
			}
			if (diagonal >= 0) values[diagonal] = -(double)sum;
			else add_csr_diagonal(row, -(double)sum, rp, ci, values);
		}
	};
	reduce_matrix(global_K_row_ptr, global_K_col_idx, global_K_values);
	reduce_matrix(global_E_row_ptr, global_E_col_idx, global_E_values);
	fixed_node.swap(reduced_fixed); fixed_value.swap(reduced_values);
	potentials.swap(reduced_potentials); global_F.swap(reduced_F);
	input_pulse_nodes.clear();
	for (auto &pulse : pulse_voltages) {
		input_pulse_nodes.push_back(pulse.nodes);
		for (int &node : pulse.nodes) node = output_node_to_dof[node];
		std::sort(pulse.nodes.begin(), pulse.nodes.end());
		pulse.nodes.erase(std::unique(pulse.nodes.begin(), pulse.nodes.end()), pulse.nodes.end());
	}
	if (circuit.enabled) {
		circuit.dof = output_node_to_dof[circuit.nodes.front()];
		for (double &v : global_K_values) v *= circuit.coordinate_scale;
		for (double &v : global_E_values) v *= circuit.coordinate_scale;
		add_csr_diagonal(circuit.dof, circuit.external_capacitance,
			global_E_row_ptr, global_E_col_idx, global_E_values);
	}
	last_dof_count = (int)nodes.size();
	printf("Equipotential reduction: %d mesh nodes -> %d degrees of freedom\n", original_count, last_dof_count);
}

void ElectrostaticAnalyzer::restore_full_system()
{
	if (output_node_to_dof.empty()) return;
	std::vector<double> full_potentials(input_nodes.size()), full_values(input_nodes.size());
	std::vector<bool> full_fixed(input_nodes.size());
	for (size_t i = 0; i < input_nodes.size(); ++i) {
		int group = output_node_to_dof[i];
		full_potentials[i] = potentials[group]; full_values[i] = fixed_value[group]; full_fixed[i] = fixed_node[group];
	}
	nodes.swap(input_nodes); input_nodes.clear();
	potentials.swap(full_potentials); fixed_value.swap(full_values); fixed_node.swap(full_fixed);
	for (size_t i = 0; i < pulse_voltages.size(); ++i) pulse_voltages[i].nodes.swap(input_pulse_nodes[i]);
	input_pulse_nodes.clear(); output_node_to_dof.clear(); circuit.dof = -1;
}

void ElectrostaticAnalyzer::record_circuit_state(double time)
{
	if (!circuit.enabled) return;
	std::vector<double> charge;
	matvec_csr(global_E_row_ptr, global_E_col_idx, global_E_values, potentials, charge);
	circuit.times.push_back(time); circuit.voltages.push_back(potentials[circuit.dof]);
	circuit.source_voltages.push_back(circuit_source_voltage(time)); circuit.charges.push_back(charge[circuit.dof]);
}

void ElectrostaticAnalyzer::write_circuit_history() const
{
	if (!circuit.enabled) return;
	size_t separator = report_output_path.find_last_of("/\\");
	std::string path = (separator == std::string::npos ? "" : report_output_path.substr(0, separator + 1)) + "circuit_history.csv";
	std::ofstream out(path);
	if (!out) electrical_error("cannot write circuit history: " + path);
	out.precision(17);
	out << "time_s,electrode_V,source_V,electrode_charge_C,branch_current_into_electrode_A\n";
	for (size_t i = 0; i < circuit.times.size(); ++i) {
		double r = circuit.times[i] < circuit.disconnect_time ? circuit.series_resistance : circuit.off_resistance;
		out << circuit.times[i] << ',' << circuit.voltages[i] << ',' << circuit.source_voltages[i] << ',' << circuit.charges[i] << ',';
		if (r > 0 && std::isfinite(r)) out << (circuit.source_voltages[i] - circuit.voltages[i]) / r;
		else if (std::isinf(r)) out << 0;
		// The ideal voltage source current is an unconstrained reaction, not V/R.
		out << '\n';
	}
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
		const size_t comment = line.find('#');
		if (comment != std::string::npos) line = line.substr(0, comment);
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
	time0 = (time_vals.size() >= 2) ? time_vals[1] : 0.0;

	if (time_vals.empty() || params.empty() || !std::isfinite(dt) || !std::isfinite(time0)) {
		fprintf(stderr, "Invalid time configuration in %s\n", filename);
		std::exit(EXIT_FAILURE);
	}
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
		} else if (bc.type == 5 && bc.value == 0.0) {
			// 零流束の自然境界は弱形式で自動的に満たされる。
			continue;
		} else {
			fprintf(stderr, "Unsupported boundary condition: type=%d value=%.6e\n", bc.type, bc.value);
			std::exit(EXIT_FAILURE);
		}
	}

	for (int i = 0; i < (int)nodes.size(); ++i) {
		if (fixed_node[i]) potentials[i] = fixed_value[i];
	}
}

void ElectrostaticAnalyzer::assemble_global_matrix()
{
	// 各四面体要素について局所行列を計算し、全体行列へ加算する。
	// 出力：K_sigma（導電率）、K_epsilon（誘電率）、F（ソース項）
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
				// Jの列は物理空間の辺ベクトルなので、勾配はJ^{-T}で変換する。
				gradN[a][r] = invJ[0][r] * dN_ref[a][0] + invJ[1][r] * dN_ref[a][1] + invJ[2][r] * dN_ref[a][2];
			}
		}

		int mat = el.material;
		if (mat < 0 || mat >= (int)materials.size()) mat = 0;
		const Material &M = materials[mat];

		double volume = std::fabs(detJ) / 6.0;
		int gid[4] = {i0, i1, i2, i3};

		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b) {
				// 異方性導電率による局所行列：
				// ke_sigma = V * (sx*Nax*Nbx + sy*Nay*Nby + sz*Naz*Nbz)
				double kxx = M.sx * gradN[a][0] * gradN[b][0];
				double kyy = M.sy * gradN[a][1] * gradN[b][1];
				double kzz = M.sz * gradN[a][2] * gradN[b][2];
				double ke_sigma = (kxx + kyy + kzz) * volume;
				global_K_rows[gid[a]].push_back(std::make_pair(gid[b], ke_sigma));

				// 時間項で使う異方性誘電率の局所行列K_epsilon。
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
	// 現在のglobal_Kとglobal_Fを使って A*phi=b を解く。
	// 導電率の桁差が大きいモデルでも安定しやすいように、規模やdtに関係なく
	// IC(0)-CG（ICCG）を使う。SSORは軽いが、強い材料コントラストでは残差が
	// 下がりにくいことがあるため、この解析では使わない。
	auto t_start = std::chrono::high_resolution_clock::now();
	const int n = (int)nodes.size();
	const bool steady_state = (dt <= 0.0);
	const bool use_ssor = false;
	const int max_iter = steady_state ? 2000 : std::max(1000, nstep * 10);
	const double tol = steady_state ? 1e-8 : 1e-10;
	use_ssor_preconditioner = use_ssor;
	last_solver_converged = false;
	last_solver_relative_residual = 1.0;

	std::vector<double> rhs = global_F;
	apply_dirichlet_constraints(rhs);

	// 対称対角スケーリングを適用する：
	//   A_s = S A S, b_s = S b, phi = S y, S_ii = 1/sqrt(A_ii)
	// CG法に必要な対称性・正定値性を維持しながら、前処理を作る前に
	// 行列係数の桁差を小さくする。過渡解析でも材料の導電率差は残るため、
	// 定常解析と同じくスケーリングしてからICCGを使う。
	std::vector<double> scale(n, 1.0);
	for (int i = 0; i < n; ++i) {
		double diag = 0.0;
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			if (global_K_col_idx[p] == i) {
				diag = global_K_values[p];
				break;
			}
		}
		if (std::isfinite(diag) && diag > 1e-30) scale[i] = 1.0 / std::sqrt(diag);
	}
	for (int i = 0; i < n; ++i) {
		for (int p = global_K_row_ptr[i]; p < global_K_row_ptr[i + 1]; ++p) {
			global_K_values[p] *= scale[i] * scale[global_K_col_idx[p]];
		}
		rhs[i] *= scale[i];
	}
	printf("Applied symmetric diagonal scaling; max_iter=%d tol=%.1e\n", max_iter, tol);

	std::vector<double> x(n, 0.0);
	for (int i = 0; i < n; ++i) x[i] = potentials[i] / scale[i];
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) x[i] = fixed_value[i] / scale[i];
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
		// IC(0)前処理を構築する。
		build_incomplete_cholesky();
		printf("Using IC(0) preconditioner for n=%d\n", n);
	}

	std::vector<double> r(n, 0.0), z(n, 0.0), p(n, 0.0), Ap(n, 0.0);
	matvec(x, Ap);
	for (int i = 0; i < n; ++i) {
		r[i] = rhs[i] - Ap[i];
		if (fixed_node[i]) r[i] = 0.0;
	}
	// 残差rに前処理を適用し、z = M^{-1}rを求める。
	apply_preconditioner(r, z);
	p = z;
	auto dot = [](const std::vector<double> &a, const std::vector<double> &b) -> double {
		double s = 0.0;
		for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
		return s;
	};

	double rnorm0 = std::sqrt(dot(r, r));
	if (rnorm0 < 1e-30) {
		potentials.resize(n);
		for (int i = 0; i < n; ++i) potentials[i] = scale[i] * x[i];
		last_solver_iterations = 0;
		last_solver_relative_residual = 0.0;
		last_solver_converged = true;
		last_solver_method = use_ssor_preconditioner ? "SSOR-CG" : "ICCG (scaled)";
		last_solver_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t_start).count();
		printf("Converged at iteration 0\n");
		return;
	}

	double rzold = dot(r, z);
	const bool log_every_cg_iteration =
		std::getenv("FEM_CG_LOG_EVERY_ITERATION") != nullptr;

	int iter_done = 0;
	for (int iter = 0; iter < max_iter; ++iter) {
		double rnorm = std::sqrt(dot(r, r));
		double rel_rnorm = rnorm / rnorm0;
		if (log_every_cg_iteration || iter < 10 || (iter + 1) % 100 == 0) {
			printf("CG iteration %d: residual=%.3e rel=%.3e\n", iter + 1, rnorm, rel_rnorm);
		}
		matvec(p, Ap);
		double denom = dot(p, Ap);
		// 正の微小値を破綻扱いすると、既に小さい残差の補正を解けない。
		if (!std::isfinite(denom) || denom <= 0.0 || !std::isfinite(rzold)) {
			printf("WARNING: CG breakdown at iteration %d\n", iter + 1);
			break;
		}

		double alpha = rzold / denom;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			if (fixed_node[i]) {
				x[i] = fixed_value[i] / scale[i];
				r[i] = 0.0;
			}
		}

		// 初期残差に対する相対残差が許容誤差未満なら収束とする。
		double rnorm2 = std::sqrt(dot(r, r));
		iter_done = iter + 1;
		last_solver_relative_residual = rnorm2 / rnorm0;
		if (last_solver_relative_residual < tol) {
			printf("Converged at iteration %d\n", iter_done);
			last_solver_converged = true;
			break;
		}

		apply_preconditioner(r, z);
		double rznew = dot(r, z);
		if (!std::isfinite(rznew) || std::fabs(rzold) < 1e-300) {
			printf("WARNING: CG preconditioner breakdown at iteration %d\n", iter + 1);
			break;
		}
		double beta = rznew / rzold;
		for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
		rzold = rznew;
	}

	potentials.resize(n);
	for (int i = 0; i < n; ++i) potentials[i] = scale[i] * x[i];
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) potentials[i] = fixed_value[i];
	}

	auto t_end = std::chrono::high_resolution_clock::now();
	double elapsed = std::chrono::duration<double>(t_end - t_start).count();
	last_solver_iterations = iter_done;
	last_solver_time = elapsed;
	last_solver_method = use_ssor_preconditioner ? "SSOR-CG" : "ICCG (scaled)";
	if (last_solver_converged) {
		printf("Solved linear system with %s: rel_residual=%.3e (%.6f s)\n",
			last_solver_method.c_str(), last_solver_relative_residual, elapsed);
	} else {
		printf("WARNING: linear solver did not converge with %s: iterations=%d rel_residual=%.3e (%.6f s)\n",
			last_solver_method.c_str(), last_solver_iterations, last_solver_relative_residual, elapsed);
		fprintf(stderr, "Analysis aborted because the linear solver did not converge.\n");
		std::exit(EXIT_FAILURE);
	}
}

void ElectrostaticAnalyzer::solve()
{
	// 支配方程式：
	//   定常解析(dt<=0)：K_sigma * phi = F
	//   過渡解析(dt>0)：(K_sigma + K_epsilon/dt) * phi^(n+1)
	//                    = F + (K_epsilon/dt) * phi^n
	const double conv_tol = 1e-8;
	printf("Starting electrostatic analysis...\n");
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;
	const bool use_pulse_voltages = !pulse_voltages.empty();
	circuit.times.clear(); circuit.voltages.clear(); circuit.source_voltages.clear(); circuit.charges.clear();
	if (circuit.enabled) {
		if (dt <= 0) electrical_error("circuit coupling requires transient analysis");
		if (std::isfinite(circuit.disconnect_time)) {
			double steps = (circuit.disconnect_time - time0) / dt;
			if (steps <= 0 || std::fabs(steps - std::round(steps)) > 1e-9)
				electrical_error("disconnect time must be after time0 and aligned with the time-step grid");
		}
	}
	for (VoltageDifferenceMonitor &monitor : voltage_difference_monitors) {
		monitor.voltage_a.clear();
		monitor.voltage_b.clear();
	}

	potentials.assign(nodes.size(), 0.0);
	if (fixed_node.size() != nodes.size()) {
		fixed_node.assign(nodes.size(), false);
		fixed_value.assign(nodes.size(), 0.0);
	}
	if (circuit.enabled) for (int node : circuit.nodes) fixed_node[node] = false;

	if (dt <= 0.0) {
		// --------------------------------------------------------------------
		// 定常解析：時間項を使わず、導電率行列による方程式を1回だけ解く。
		// --------------------------------------------------------------------
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
		prepare_reduced_system();
		solve_linear_system();
		for (VoltageDifferenceMonitor &monitor : voltage_difference_monitors) {
			monitor.voltage_a.push_back(potentials[monitor_dof(monitor.node_a)]);
			monitor.voltage_b.push_back(potentials[monitor_dof(monitor.node_b)]);
		}
		final_step = 1;
		final_time = 0.0;
		final_max_diff = 0.0;
		step_solver_iterations.clear();
		step_solver_times.clear();
		step_solver_iterations.push_back(last_solver_iterations);
		step_solver_times.push_back(last_solver_time);
		write_potential_distribution((step_output_prefix + "001").c_str());
		restore_full_system();
		printf("Analysis completed.\n");
		write_analysis_report(report_output_path.c_str());
		return;
	}

	// ------------------------------------------------------------------------
	// 過渡解析の準備：境界条件、行列組立、RCM並べ替え
	// ------------------------------------------------------------------------
	apply_boundary_conditions();
	if (use_pulse_voltages) {
		apply_pulse_voltages(time0);
	}
	apply_circuit_boundary(time0);
	if (circuit.enabled) {
		bool grounded = false;
		for (size_t i = 0; i < fixed_node.size(); ++i)
			if (fixed_node[i] && fixed_value[i] == 0.0 &&
				!std::binary_search(circuit.nodes.begin(), circuit.nodes.end(), (int)i)) grounded = true;
		if (!grounded) electrical_error("circuit model needs an explicit 0 V reference outside the driven electrode");
	}
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
	prepare_reduced_system();
	const int n = (int)nodes.size();
	if (n > 0) {
		// t=0で既にHighのパルスも、自由節点に初期電荷を作らないようにする。
		// 境界電位に整合する誘電場 K_epsilon * phi(0) = 0 を自由節点で解く。
		// 導電定常場ではなく、瞬時印加直後の無自由電荷状態を初期値とする。
		std::vector<int> saved_row_ptr = global_K_row_ptr;
		std::vector<int> saved_col_idx = global_K_col_idx;
		std::vector<double> saved_values = global_K_values;
		std::vector<double> saved_F = global_F;
		if (circuit.enabled && circuit.steady_initial_state) {
			add_circuit_terms(time0);
			printf("Initializing conductive steady state at t=%.6e s.\n", time0);
		} else {
			global_K_row_ptr = global_E_row_ptr;
			global_K_col_idx = global_E_col_idx;
			global_K_values = global_E_values;
			global_F.assign(n, 0.0);
			printf("Initializing charge-free dielectric field at t=%.6e s.\n", time0);
		}
		solve_linear_system();
		global_K_row_ptr.swap(saved_row_ptr);
		global_K_col_idx.swap(saved_col_idx);
		global_K_values.swap(saved_values);
		global_F.swap(saved_F);
		if (!last_solver_converged) {
			fprintf(stderr, "Failed to initialize dielectric field; check permittivity and boundary constraints.\n");
			std::exit(EXIT_FAILURE);
		}
	}
	record_circuit_state(time0);
	std::vector<int> perm(n, 0);
	std::vector<int> inv_perm(n, 0);
	if (n > 0) {
		// RCM（Reverse Cuthill-McKee）法で節点を並べ替える。
		// 非ゼロ成分を対角線付近へ集め、行列の帯幅と前処理コストを抑える。
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

	for (PulseVoltageEntry &pulse : pulse_voltages) {
		for (int &node : pulse.nodes) {
			if (0 <= node && node < n) node = inv_perm[node];
		}
	}
	if (circuit.enabled) circuit.dof = inv_perm[circuit.dof];

	std::vector<double> base_F = global_F;
	// RCM後のK_sigmaを「過渡解析の元行列」として保存する。
	// solve_linear_system()はDirichlet条件やスケーリングをglobal_Kへ直接埋め込むため、
	// 毎ステップここから復元してから A = K_sigma + K_epsilon/dt を作る。
	// これにより、K_epsilon/dt がステップごとに二重・三重に加算されることを防ぐ。
	const std::vector<int> base_K_row_ptr = global_K_row_ptr;
	const std::vector<int> base_K_col_idx = global_K_col_idx;
	const std::vector<double> base_K_values = global_K_values;

	auto write_original_order_distribution = [&](const char *filename) {
		// 計算中のpotentialsはRCM後の節点順なので、tempa.datへ出す直前だけ
		// 入力メッシュの節点順へ戻す。内部の計算状態はRCM順のまま維持する。
		std::vector<double> permuted_potentials = potentials;
		std::vector<double> original_order_potentials(n, 0.0);
		for (int i = 0; i < n; ++i) {
			original_order_potentials[perm[i]] = permuted_potentials[i];
		}
		potentials.swap(original_order_potentials);
		write_potential_distribution(filename);
		potentials.swap(permuted_potentials);
	};

	printf("Transient solve: dt=%.3e, nstep=%d\n", dt, nstep);
	if (use_pulse_voltages) {
		printf("Pulse voltage enabled: early convergence disabled; running all %d steps.\n", nstep);
	}
	step_solver_iterations.clear();
	step_solver_times.clear();

	std::vector<double> potentials_prev = potentials;
	std::vector<int> step_row_ptr, step_col_idx;
	std::vector<double> step_values, eps_phi, potentials_old;
	for (int step = 0; step < nstep; ++step) {
		// --------------------------------------------------------------------
		// Backward Eulerによる時間ステップ
		//   A   = K_sigma + K_epsilon/dt
		//   rhs = F + (K_epsilon/dt) * phi^n
		// --------------------------------------------------------------------
		potentials_old = potentials_prev;
		if (use_pulse_voltages) {
			apply_pulse_voltages(time0 + (step + 1) * dt);
		}
		double network_time = circuit_step_time(time0 + (step + 1) * dt);
		apply_circuit_boundary(network_time);
		global_K_row_ptr = base_K_row_ptr;
		global_K_col_idx = base_K_col_idx;
		global_K_values = base_K_values;
		global_F = base_F;
		add_circuit_terms(network_time);
		build_time_step_matrix(1.0 / dt, step_row_ptr, step_col_idx, step_values);

		global_K_row_ptr.swap(step_row_ptr);
		global_K_col_idx.swap(step_col_idx);
		global_K_values.swap(step_values);

		matvec_csr(global_E_row_ptr, global_E_col_idx, global_E_values, potentials_prev, eps_phi);
		for (size_t i = 0; i < global_F.size(); ++i) {
			global_F[i] += eps_phi[i] / dt;
		}

		potentials = potentials_prev;
		if (use_pulse_voltages) {
			apply_pulse_voltages(time0 + (step + 1) * dt);
		}
		apply_circuit_boundary(network_time);
		solve_linear_system();
		potentials_prev = potentials;
		record_circuit_state(time0 + (step + 1) * dt);
		for (VoltageDifferenceMonitor &monitor : voltage_difference_monitors) {
			monitor.voltage_a.push_back(potentials[inv_perm[monitor_dof(monitor.node_a)]]);
			monitor.voltage_b.push_back(potentials[inv_perm[monitor_dof(monitor.node_b)]]);
		}
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
			step + 1, nstep, time0 + (step + 1) * dt, vmin, vmax, max_diff);

		char step_number[32];
		std::snprintf(step_number, sizeof(step_number), "%03d", step + 1);
		const std::string outname = step_output_prefix + step_number;
		write_original_order_distribution(outname.c_str());

		final_step = step + 1;
		final_time = time0 + (step + 1) * dt;
		final_max_diff = max_diff;

		if (!use_pulse_voltages && !circuit.enabled && max_diff < conv_tol) {
			printf("Converged at step %d: max_diff=%.6e < tol=%.6e\n", step + 1, max_diff, conv_tol);
			converged_early = true;
			break;
		}
	}

	// RCMで変更した節点順を、入力メッシュの節点順へ戻す。
	std::vector<double> original_potentials(n, 0.0);
	for (int i = 0; i < n; ++i) {
		original_potentials[perm[i]] = potentials[i];
	}
	potentials.swap(original_potentials);
	// 再度solve()を呼べるように、境界条件の節点番号も入力順へ戻す。
	std::vector<bool> original_fixed_node(n, false);
	std::vector<double> original_fixed_value(n, 0.0);
	for (int i = 0; i < n; ++i) {
		original_fixed_node[perm[i]] = fixed_node[i];
		original_fixed_value[perm[i]] = fixed_value[i];
	}
	fixed_node.swap(original_fixed_node);
	fixed_value.swap(original_fixed_value);
	for (PulseVoltageEntry &pulse : pulse_voltages) {
		for (int &node : pulse.nodes) node = perm[node];
	}
	if (circuit.enabled) circuit.dof = perm[circuit.dof];
	restore_full_system();
	write_circuit_history();

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

	fprintf(fp, "Electrostatic FEM Analysis Summary / 静電場FEM解析結果\n");
	fprintf(fp, "=====================================================\n\n");
	fprintf(fp, "このファイルの見方\n");
	fprintf(fp, "------------------\n");
	fprintf(fp, "- dt <= 0 は定常解析、dt > 0 はBackward Euler法による過渡解析です。\n");
	fprintf(fp, "- last_solver_converged=yes は、線形ソルバーの相対残差が許容値未満になったことを表します。\n");
	fprintf(fp, "- converged_early=yes は、時間ステップ間の最大電位差が収束許容値未満になったことを表します。\n");
	fprintf(fp, "- 電位の単位はV、時間の単位はs、solve_timeの単位はsです。\n\n");
	fprintf(fp, "Mesh\n");
	fprintf(fp, "----\n");
	fprintf(fp, "解析に使用した四面体メッシュの規模です。\n");
	fprintf(fp, "Nodes: %d\n", (int)nodes.size());
	fprintf(fp, "Elements: %d\n\n", (int)elements.size());
	fprintf(fp, "axial_symmetry: %s\n", axial_symmetry ? "z (equal potential along z)" : "disabled");
	fprintf(fp, "solver_degrees_of_freedom: %d\n", last_dof_count);
	if (circuit.enabled) {
		fprintf(fp, "circuit_electrode_material: %d\n", circuit.material_id + 1);
		fprintf(fp, "circuit_series_resistance_ohm: %.12e\n", circuit.series_resistance);
		fprintf(fp, "circuit_off_resistance_ohm: %.12e\n", circuit.off_resistance);
		fprintf(fp, "circuit_disconnect_time_s: %.12e\n", circuit.disconnect_time);
		fprintf(fp, "circuit_external_capacitance_F: %.12e\n", circuit.external_capacitance);
		fprintf(fp, "circuit_coordinate_scale_to_m: %.12e\n", circuit.coordinate_scale);
		fprintf(fp, "initial_state: %s\n", circuit.steady_initial_state ? "conductive steady state" : "charge-free dielectric field");
		fprintf(fp, "circuit_trace: circuit_history.csv (includes initial state)\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "Analysis Conditions\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "時間積分条件と解析全体の終了状態です。\n");
	fprintf(fp, "analysis_type: %s\n", dt <= 0.0 ? "steady-state（定常解析）" : "transient（過渡解析）");
	fprintf(fp, "dt: %.6e\n", dt);
	fprintf(fp, "time0: %.6e\n", time0);
	fprintf(fp, "nstep(max): %d\n", nstep);
	fprintf(fp, "convergence_tolerance: %.6e\n", 1e-8);
	fprintf(fp, "final_step: %d\n", final_step);
	fprintf(fp, "final_time: %.6e\n", final_time);
	fprintf(fp, "final_max_diff: %.6e\n", final_max_diff);
	fprintf(fp, "converged_early: %s\n\n", converged_early ? "yes" : "no");
	fprintf(fp, "solver_method: %s\n", last_solver_method.c_str());
	fprintf(fp, "  ICCG: 不完全コレスキー前処理付き共役勾配法\n");
	fprintf(fp, "  scaled: 行列係数の桁差を抑える対称対角スケーリングを適用\n");
	fprintf(fp, "last_solver_converged: %s\n", last_solver_converged ? "yes" : "no");
	fprintf(fp, "last_solver_relative_residual: %.6e\n", last_solver_relative_residual);
	fprintf(fp, "last_solver_iterations: %d\n", last_solver_iterations);
	fprintf(fp, "last_solver_time_sec: %.6e\n\n", last_solver_time);

	fprintf(fp, "Solver Iterations Per Step\n");
	fprintf(fp, "-------------------------\n");
	fprintf(fp, "各時間ステップにおける線形ソルバーの反復回数と求解時間です。\n");
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

	fprintf(fp, "Voltage Differences Per Step\n");
	fprintf(fp, "----------------------------\n");
	fprintf(fp, "voltage_difference.datで指定した2接点間の電位差です。difference_V = V(node_b) - V(node_a)です。\n");
	if (voltage_difference_monitors.empty()) {
		fprintf(fp, "(no voltage difference monitors)\n\n");
	} else {
		for (const VoltageDifferenceMonitor &monitor : voltage_difference_monitors) {
			fprintf(fp, "label: %s\n", monitor.label.c_str());
			fprintf(fp, "node_a: %d\n", monitor.node_a + 1);
			fprintf(fp, "node_b: %d\n", monitor.node_b + 1);
			fprintf(fp, "step,time_s,node_a_V,node_b_V,difference_V,absolute_difference_V\n");
			size_t count = std::min(monitor.voltage_a.size(), monitor.voltage_b.size());
			for (size_t i = 0; i < count; ++i) {
				double difference = monitor.voltage_b[i] - monitor.voltage_a[i];
				double time = dt > 0.0 ? time0 + (i + 1) * dt : 0.0;
				fprintf(fp, "%zu,%.6e,%.6e,%.6e,%.6e,%.6e\n",
					i + 1, time, monitor.voltage_a[i], monitor.voltage_b[i],
					difference, std::fabs(difference));
			}
			fprintf(fp, "\n");
		}
	}

	fprintf(fp, "Boundary Conditions\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "sin.datで節点番号を指定した境界条件です。type 1/4はnode1とnode2を固定し、type 5の値0は零流束の自然境界です。その他の条件は未対応のためエラーになります。\n");
	fprintf(fp, "count: %d\n", (int)boundaries.size());
	for (size_t i = 0; i < boundaries.size(); ++i) {
		const BoundaryCondition &bc = boundaries[i];
		fprintf(fp, "%zu: node1=%d node2=%d type=%d value=%.6e\n",
			i + 1, bc.node1 + 1, bc.node2 + 1, bc.type, bc.value);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Fixed Materials (tmate.dat)\n");
	fprintf(fp, "---------------------------\n");
	fprintf(fp, "指定材料に属する節点へ適用した固定電位です。valueの単位はVです。\n");
	fprintf(fp, "count: %d\n", (int)fixed_materials.size());
	for (size_t i = 0; i < fixed_materials.size(); ++i) {
		const FixedMaterialEntry &fm = fixed_materials[i];
		fprintf(fp, "%zu: material_id=%d value=%.6e\n", i + 1, fm.material_id + 1, fm.value);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Pulse Voltages (pulse.dat)\n");
	fprintf(fp, "--------------------------\n");
	fprintf(fp, "台形または矩形パルス電位です。回路対象材料では電源側の波形として使い、電極電位は回路とFEMから計算します。その他の材料では固定電位として適用します。早期収束判定は無効です。\n");
	fprintf(fp, "format: material_id v_low v_high t_rise t_high t_fall t_low [phase]\n");
	fprintf(fp, "count: %d\n", (int)pulse_voltages.size());
	for (size_t i = 0; i < pulse_voltages.size(); ++i) {
		const PulseVoltageEntry &pulse = pulse_voltages[i];
		const double period = pulse.t_rise + pulse.t_high + pulse.t_fall + pulse.t_low;
		const bool is_square = pulse.t_rise == 0.0 && pulse.t_fall == 0.0;
		const double high_start = pulse.t_rise;
		const double high_end = high_start + pulse.t_high;
		const double fall_end = high_end + pulse.t_fall;
		const double high_ratio = period > 0.0 ? 100.0 * pulse.t_high / period : 0.0;
		fprintf(fp, "%zu: material_id=%d v_low=%.6e v_high=%.6e t_rise=%.6e t_high=%.6e t_fall=%.6e t_low=%.6e phase=%.6e nodes=%d\n",
			i + 1, pulse.material_id + 1, pulse.v_low, pulse.v_high, pulse.t_rise,
			pulse.t_high, pulse.t_fall, pulse.t_low, pulse.phase, (int)pulse.nodes.size());
		fprintf(fp, "   waveform_type: %s\n", is_square ? "square（矩形波）" : "trapezoid（台形波）");
		fprintf(fp, "   period: %.6e s (%.6f ns)\n", period, period * 1.0e9);
		if (period > 0.0) {
			fprintf(fp, "   frequency: %.6e Hz\n", 1.0 / period);
		}
		fprintf(fp, "   %s: %.3f %%\n",
			is_square ? "duty_cycle" : "high_hold_ratio", high_ratio);
		fprintf(fp, "   cycle_timeline: rise=[0, %.6f) ns, high=[%.6f, %.6f) ns, fall=[%.6f, %.6f) ns, low=[%.6f, %.6f) ns\n",
			high_start * 1.0e9, high_start * 1.0e9, high_end * 1.0e9,
			high_end * 1.0e9, fall_end * 1.0e9, fall_end * 1.0e9, period * 1.0e9);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Material Properties\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "材料ごとの異方性物性です。sigma=(x,y,z)は導電率、epsilon=(x,y,z)は誘電率です。\n");
	fprintf(fp, "count: %d\n", (int)materials.size());
	for (size_t i = 0; i < materials.size(); ++i) {
		const Material &m = materials[i];
		fprintf(fp, "%zu: sigma=(%.6e, %.6e, %.6e) epsilon=(%.6e, %.6e, %.6e)\n",
			i + 1, m.sx, m.sy, m.sz, m.ex, m.ey, m.ez);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Material Block Ranges\n");
	fprintf(fp, "---------------------\n");
	fprintf(fp, "要素番号startからendまでに割り当てられた材料IDです（番号は1始まり、両端を含む）。\n");
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

		// 1つの値を12文字以内に収める。
		// 指数部が3桁必要な値は、符号付きゼロとして出力する。
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
	size_t output_count = output_node_to_dof.empty() ? potentials.size() : output_node_to_dof.size();
	for (size_t i = 0; i < output_count; ++i) {
		double value = potentials[output_node_to_dof.empty() ? i : (size_t)output_node_to_dof[i]];
		if (!std::isfinite(value)) electrical_error("nonfinite nodal potential during output");
		fprintf(fp, "%s", format_e12_4(value).c_str());
		col++;
		if (col == 6 || i + 1 == output_count) {
			fprintf(fp, "\n");
			col = 0;
		}
	}

	fclose(fp);
	printf("Wrote %s: %d nodal potentials\n", filename, (int)output_count);
}
