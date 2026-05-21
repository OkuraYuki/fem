#include "electrostatic_analyzer.h"

#include <cstdio>
#include <cctype>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

ElectrostaticAnalyzer::ElectrostaticAnalyzer()
	: dt(0.0), nstep(200), unit(1.0)
{
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;
	last_cg_iterations = 0;
	last_cg_initial_residual = 0.0;
	last_cg_final_residual = 0.0;
	last_cg_converged = false;
	step_output_prefix = "tempa.dat";
	report_output_path = "analysis_summary.txt";
	fixed_node_count = 0;
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
	const int n = (int)nodes.size();
	std::vector<std::vector<std::pair<int, double>>> rows(n);
	for (int r = 0; r < n; ++r) {
		rows[r] = global_K_rows[r];
		for (const std::pair<int, double> &entry : global_E_rows[r]) {
			rows[r].push_back(std::make_pair(entry.first, entry.second * scale_eps));
		}
	}
	build_sparse_from_rows(rows, row_ptr, col_idx, values);
}

void ElectrostaticAnalyzer::build_sparse_matrix()
{
	build_sparse_from_rows(global_K_rows, global_K_row_ptr, global_K_col_idx, global_K_values);
	build_sparse_from_rows(global_E_rows, global_E_row_ptr, global_E_col_idx, global_E_values);
}

void ElectrostaticAnalyzer::apply_dirichlet_constraints(std::vector<double> &rhs)
{
	const int n = (int)nodes.size();
	for (int i = 0; i < n; ++i) {
		if (!fixed_node[i]) continue;

		for (int r = 0; r < n; ++r) {
			for (int p = global_K_row_ptr[r]; p < global_K_row_ptr[r + 1]; ++p) {
				if (global_K_col_idx[p] != i) continue;
				if (r != i) {
					rhs[r] -= global_K_values[p] * fixed_value[i];
				}
				global_K_values[p] = 0.0;
			}
		}

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
			global_K_col_idx.insert(global_K_col_idx.begin() + global_K_row_ptr[i + 1], i);
			global_K_values.insert(global_K_values.begin() + global_K_row_ptr[i + 1], 1.0);
			for (int r = i + 1; r <= n; ++r) {
				++global_K_row_ptr[r];
			}
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

void ElectrostaticAnalyzer::collect_input_quality_stats()
{
	input_quality.invalid_connectivity_elements = 0;
	input_quality.invalid_connectivity_references = 0;
	input_quality.isolated_nodes = 0;
	input_quality.invalid_boundary_nodes = 0;
	input_quality.unassigned_material_elements = 0;

	std::vector<int> node_usage(nodes.size(), 0);
	for (size_t e = 0; e < elements.size(); ++e) {
		bool valid = true;
		for (int k = 0; k < 4; ++k) {
			const int n = elements[e].nodes[k];
			if (n < 0 || n >= (int)nodes.size()) {
				++input_quality.invalid_connectivity_references;
				valid = false;
			} else {
				++node_usage[n];
			}
		}
		if (!valid) ++input_quality.invalid_connectivity_elements;
		if (elements[e].material < 0 || elements[e].material >= (int)materials.size()) {
			++input_quality.unassigned_material_elements;
		}
	}

	for (size_t i = 0; i < node_usage.size(); ++i) {
		if (node_usage[i] == 0) ++input_quality.isolated_nodes;
	}

	for (size_t i = 0; i < boundaries.size(); ++i) {
		const BoundaryCondition &bc = boundaries[i];
		if (bc.node1 < 0 || bc.node1 >= (int)nodes.size()) ++input_quality.invalid_boundary_nodes;
		if (bc.node2 < 0 || bc.node2 >= (int)nodes.size()) ++input_quality.invalid_boundary_nodes;
	}
}

void ElectrostaticAnalyzer::collect_pre_solve_summary()
{
	fixed_node_count = 0;
	for (size_t i = 0; i < fixed_node.size(); ++i) {
		if (fixed_node[i]) ++fixed_node_count;
	}

	std::map<int, int> bc_hist;
	for (size_t i = 0; i < boundaries.size(); ++i) {
		++bc_hist[boundaries[i].type];
	}
	boundary_type_counts.clear();
	for (std::map<int, int>::const_iterator it = bc_hist.begin(); it != bc_hist.end(); ++it) {
		boundary_type_counts.push_back(std::make_pair(it->first, it->second));
	}

	material_element_counts.assign(materials.size(), 0);
	for (size_t e = 0; e < elements.size(); ++e) {
		int mat = elements[e].material;
		if (mat >= 0 && mat < (int)material_element_counts.size()) ++material_element_counts[mat];
	}

	printf("Pre-solve summary: fixed_nodes=%d, boundary_count=%d, material_count=%d\n",
		fixed_node_count, (int)boundaries.size(), (int)materials.size());
}

MatrixDiagnostics ElectrostaticAnalyzer::analyze_matrix(const std::vector<int> &row_ptr,
	const std::vector<int> &col_idx, const std::vector<double> &values) const
{
	MatrixDiagnostics stats;
	stats.nnz = (int)values.size();
	stats.min_abs_diag = std::numeric_limits<double>::infinity();

	const int n = (int)nodes.size();
	for (int r = 0; r < n; ++r) {
		const int row_nnz = row_ptr[r + 1] - row_ptr[r];
		if (row_nnz > stats.max_row_nnz) stats.max_row_nnz = row_nnz;

		bool has_diag = false;
		double diag_abs = 0.0;
		for (int p = row_ptr[r]; p < row_ptr[r + 1]; ++p) {
			double av = std::fabs(values[p]);
			if (av > stats.max_abs_value) stats.max_abs_value = av;
			if (col_idx[p] == r) {
				has_diag = true;
				diag_abs = av;
			}
		}

		if (!has_diag || diag_abs <= 0.0) {
			++stats.zero_diag_count;
		} else {
			if (diag_abs < stats.min_abs_diag) stats.min_abs_diag = diag_abs;
			if (diag_abs > stats.max_abs_diag) stats.max_abs_diag = diag_abs;
		}
	}

	if (!std::isfinite(stats.min_abs_diag)) stats.min_abs_diag = 0.0;
	return stats;
}

void ElectrostaticAnalyzer::log_matrix_diagnostics(const char *label, int step) const
{
	printf("[%s][step=%d] nnz=%d, max_row_nnz=%d, zero_diag=%d, min|diag|=%.6e, max|diag|=%.6e, max|A|=%.6e\n",
		label, step, last_matrix_diagnostics.nnz, last_matrix_diagnostics.max_row_nnz,
		last_matrix_diagnostics.zero_diag_count, last_matrix_diagnostics.min_abs_diag,
		last_matrix_diagnostics.max_abs_diag, last_matrix_diagnostics.max_abs_value);
}

void ElectrostaticAnalyzer::read_mesh(const char *filename)
{
	input_quality.invalid_connectivity_elements = 0;
	input_quality.invalid_connectivity_references = 0;
	input_quality.isolated_nodes = 0;
	input_quality.negative_jacobian_elements = 0;
	input_quality.degenerate_elements = 0;
	input_quality.material_fallback_elements = 0;

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

	elements.assign(num_elements, Element{{0, 0, 0, 0}, -1});
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
	printf("Read mesh: %d nodes, %d elements\n", num_nodes, num_elements);
}

void ElectrostaticAnalyzer::read_material_and_bc(const char *filename)
{
	input_quality.invalid_boundary_nodes = 0;
	input_quality.invalid_material_block_ranges = 0;
	input_quality.invalid_material_ids = 0;
	input_quality.unassigned_material_elements = 0;

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
		if (start > end || start < 1 || end > (int)elements.size()) {
			++input_quality.invalid_material_block_ranges;
		}
		if (mat_id < 0 || mat_id >= ntmate) {
			++input_quality.invalid_material_ids;
		}
		material_ranges.push_back(MaterialBlockRange{start, end, mat_id});
		for (int e = start - 1; e <= end - 1; ++e) {
			if (0 <= e && e < (int)elements.size()) {
				elements[e].material = mat_id;
			}
		}
	}

	if (ntmblk == 0 && ntmate > 0) {
		for (size_t e = 0; e < elements.size(); ++e) {
			elements[e].material = 0;
		}
		printf("No material block ranges found; fallback to material_id=1 for all elements.\n");
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
	fixed_materials.clear();

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
	const int n = (int)nodes.size();
	global_K_rows.assign(n, {});
	global_E_rows.assign(n, {});
	global_F.assign(n, 0.0);
	input_quality.negative_jacobian_elements = 0;
	input_quality.degenerate_elements = 0;
	input_quality.material_fallback_elements = 0;

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

		if (detJ < 0.0) ++input_quality.negative_jacobian_elements;
		if (std::fabs(detJ) < 1e-20) {
			++input_quality.degenerate_elements;
			continue;
		}

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
		if (mat < 0 || mat >= (int)materials.size()) {
			mat = 0;
			++input_quality.material_fallback_elements;
		}
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

	printf("Assembled sparse global matrix: %dx%d, nnz=%d\n", n, n, (int)global_K_values.size());
}

void ElectrostaticAnalyzer::solve_linear_system()
{
	const int n = (int)nodes.size();
	const int max_iter = std::max(50, nstep);
	const double tol = 1e-10;

	std::vector<double> rhs = global_F;
	apply_dirichlet_constraints(rhs);

	std::vector<double> x = potentials;
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) x[i] = fixed_value[i];
	}

	std::vector<double> r(n, 0.0), p(n, 0.0), Ap(n, 0.0);
	matvec(x, Ap);
	for (int i = 0; i < n; ++i) {
		r[i] = rhs[i] - Ap[i];
		if (fixed_node[i]) r[i] = 0.0;
		p[i] = r[i];
	}

	auto dot = [](const std::vector<double> &a, const std::vector<double> &b) -> double {
		double s = 0.0;
		for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
		return s;
	};

	double rsold = dot(r, r);
	last_cg_initial_residual = std::sqrt(rsold);
	last_cg_final_residual = last_cg_initial_residual;
	last_cg_iterations = 0;
	last_cg_converged = false;
	printf("CG start: initial_residual=%.6e, tol=%.1e, max_iter=%d\n", last_cg_initial_residual, tol, max_iter);
	if (rsold < tol * tol) {
		potentials = x;
		last_cg_converged = true;
		printf("Converged at iteration 0\n");
		printf("CG summary: iterations=0, final_residual=%.6e, converged=yes\n", last_cg_final_residual);
		return;
	}

	int iter_done = 0;
	for (int iter = 0; iter < max_iter; ++iter) {
		matvec(p, Ap);
		double denom = dot(p, Ap);
		if (std::fabs(denom) < 1e-30) break;

		double alpha = rsold / denom;
		for (int i = 0; i < n; ++i) {
			x[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
			if (fixed_node[i]) {
				x[i] = fixed_value[i];
				r[i] = 0.0;
			}
		}

		double rsnew = dot(r, r);
		iter_done = iter + 1;
		last_cg_final_residual = std::sqrt(rsnew);
		if (iter_done == 1 || iter_done % 50 == 0) {
			printf("CG iter %d: residual=%.6e\n", iter_done, last_cg_final_residual);
		}
		if (std::sqrt(rsnew) < tol) {
			last_cg_converged = true;
			printf("Converged at iteration %d\n", iter_done);
			break;
		}

		double beta = rsnew / rsold;
		for (int i = 0; i < n; ++i) {
			p[i] = r[i] + beta * p[i];
		}
		rsold = rsnew;
	}

	potentials = x;
	last_cg_iterations = iter_done;
	for (int i = 0; i < n; ++i) {
		if (fixed_node[i]) potentials[i] = fixed_value[i];
	}

	printf("CG summary: iterations=%d, final_residual=%.6e, converged=%s\n",
		last_cg_iterations, last_cg_final_residual, last_cg_converged ? "yes" : "no");
}

void ElectrostaticAnalyzer::solve()
{
	const double conv_tol = 1e-8;
	printf("Starting electrostatic analysis...\n");
	final_step = 0;
	final_time = 0.0;
	final_max_diff = 0.0;
	converged_early = false;
	boundary_type_counts.clear();
	material_element_counts.clear();
	fixed_node_count = 0;
	last_matrix_diagnostics = MatrixDiagnostics();

	potentials.assign(nodes.size(), 0.0);
	if (fixed_node.size() != nodes.size()) {
		fixed_node.assign(nodes.size(), false);
		fixed_value.assign(nodes.size(), 0.0);
	}

	if (dt <= 0.0) {
		printf("Invalid time step: dt must be positive for transient analysis.\n");
		return;
	}

	collect_input_quality_stats();
	printf("Input quality: invalid_connectivity_elements=%d, invalid_connectivity_refs=%d, isolated_nodes=%d, unassigned_material_elements=%d, invalid_bc_nodes=%d, invalid_material_ranges=%d, invalid_material_ids=%d\n",
		input_quality.invalid_connectivity_elements,
		input_quality.invalid_connectivity_references,
		input_quality.isolated_nodes,
		input_quality.unassigned_material_elements,
		input_quality.invalid_boundary_nodes,
		input_quality.invalid_material_block_ranges,
		input_quality.invalid_material_ids);

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
	collect_pre_solve_summary();
	assemble_global_matrix();
	std::vector<double> base_F = global_F;
	printf("Transient solve: dt=%.3e, nstep=%d\n", dt, nstep);

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
		if (step == 0 || step + 1 == nstep) {
			last_matrix_diagnostics = analyze_matrix(global_K_row_ptr, global_K_col_idx, global_K_values);
			log_matrix_diagnostics("step_matrix", step + 1);
		}

		global_F = base_F;
		std::vector<double> eps_phi;
		matvec_csr(global_E_row_ptr, global_E_col_idx, global_E_values, potentials_prev, eps_phi);
		for (size_t i = 0; i < global_F.size(); ++i) {
			global_F[i] += eps_phi[i] / dt;
		}

		potentials = potentials_prev;
		solve_linear_system();
		potentials_prev = potentials;

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

	fprintf(fp, "Input Quality Checks\n");
	fprintf(fp, "--------------------\n");
	fprintf(fp, "invalid_connectivity_elements: %d\n", input_quality.invalid_connectivity_elements);
	fprintf(fp, "invalid_connectivity_references: %d\n", input_quality.invalid_connectivity_references);
	fprintf(fp, "isolated_nodes: %d\n", input_quality.isolated_nodes);
	fprintf(fp, "invalid_boundary_nodes: %d\n", input_quality.invalid_boundary_nodes);
	fprintf(fp, "invalid_material_block_ranges: %d\n", input_quality.invalid_material_block_ranges);
	fprintf(fp, "invalid_material_ids: %d\n", input_quality.invalid_material_ids);
	fprintf(fp, "unassigned_material_elements: %d\n", input_quality.unassigned_material_elements);
	fprintf(fp, "negative_jacobian_elements: %d\n", input_quality.negative_jacobian_elements);
	fprintf(fp, "degenerate_elements: %d\n", input_quality.degenerate_elements);
	fprintf(fp, "material_fallback_elements: %d\n\n", input_quality.material_fallback_elements);

	fprintf(fp, "Analysis Conditions\n");
	fprintf(fp, "-------------------\n");
	fprintf(fp, "dt: %.6e\n", dt);
	fprintf(fp, "nstep(max): %d\n", nstep);
	fprintf(fp, "convergence_tolerance: %.6e\n", 1e-8);
	fprintf(fp, "final_step: %d\n", final_step);
	fprintf(fp, "final_time: %.6e\n", final_time);
	fprintf(fp, "final_max_diff: %.6e\n", final_max_diff);
	fprintf(fp, "converged_early: %s\n\n", converged_early ? "yes" : "no");

	fprintf(fp, "Pre-solve Summary\n");
	fprintf(fp, "-----------------\n");
	fprintf(fp, "fixed_nodes: %d\n", fixed_node_count);
	fprintf(fp, "boundary_type_count: %d\n", (int)boundary_type_counts.size());
	for (size_t i = 0; i < boundary_type_counts.size(); ++i) {
		fprintf(fp, "boundary_type_%d: %d\n", boundary_type_counts[i].first, boundary_type_counts[i].second);
	}
	fprintf(fp, "materials_with_elements: %d\n", (int)material_element_counts.size());
	for (size_t i = 0; i < material_element_counts.size(); ++i) {
		fprintf(fp, "material_%zu_elements: %d\n", i + 1, material_element_counts[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Linear Solver Diagnostics\n");
	fprintf(fp, "-------------------------\n");
	fprintf(fp, "cg_iterations_last_step: %d\n", last_cg_iterations);
	fprintf(fp, "cg_initial_residual_last_step: %.6e\n", last_cg_initial_residual);
	fprintf(fp, "cg_final_residual_last_step: %.6e\n", last_cg_final_residual);
	fprintf(fp, "cg_converged_last_step: %s\n\n", last_cg_converged ? "yes" : "no");

	fprintf(fp, "Matrix Diagnostics (latest)\n");
	fprintf(fp, "---------------------------\n");
	fprintf(fp, "nnz: %d\n", last_matrix_diagnostics.nnz);
	fprintf(fp, "max_row_nnz: %d\n", last_matrix_diagnostics.max_row_nnz);
	fprintf(fp, "zero_diag_count: %d\n", last_matrix_diagnostics.zero_diag_count);
	fprintf(fp, "min_abs_diag: %.6e\n", last_matrix_diagnostics.min_abs_diag);
	fprintf(fp, "max_abs_diag: %.6e\n", last_matrix_diagnostics.max_abs_diag);
	fprintf(fp, "max_abs_value: %.6e\n\n", last_matrix_diagnostics.max_abs_value);

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
		char temp[64];
		// Output with 4 decimal places in scientific notation
		snprintf(buf, sizeof(buf), "%.4E", value);
		
		// Find and adjust exponent to 2 digits (Fortran E12.4 format)
		char *e_pos = strchr(buf, 'E');
		if (e_pos) {
			char exp_sign = *(e_pos + 1);
			int exp_val = std::abs(atoi(e_pos + 2));
			// Limit exponent to 2 digits for Fortran compatibility
			if (exp_val > 99) exp_val = exp_val % 100;
			snprintf(temp, sizeof(temp), "%.*sE%c%02d", (int)(e_pos - buf), buf, exp_sign, exp_val);
			snprintf(buf, sizeof(buf), "%s", temp);
		}
		
		// Right-align to 12 characters
		std::string result = buf;
		while ((int)result.length() < 12) {
			result = " " + result;
		}
		// Truncate to 12 if longer
		if ((int)result.length() > 12) {
			result = result.substr(0, 12);
		}
		return result;
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
