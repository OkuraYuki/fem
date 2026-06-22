#ifndef ELECTROSTATIC_ANALYZER_H
#define ELECTROSTATIC_ANALYZER_H

#include <vector>
#include <utility>
#include <string>

struct Node {
	double x, y, z;
};

struct Element {
	int nodes[4];      // 4節点四面体
	int material;      // 0-based material id
};

struct Material {
	double sx, sy, sz; // 導電率
	double ex, ey, ez; // 誘電率
	double ro, cv;     // 互換用ダミー
};

struct BoundaryCondition {
	int node1, node2;  // 0-based node id
	int type;          // 1/4: Dirichlet, 5: Neumann(自然境界)
	double value;
};

struct MaterialBlockRange {
	int start_element;
	int end_element;
	int material_id;
};

struct FixedMaterialEntry {
	int material_id;
	double value;
};

class ElectrostaticAnalyzer {
private:
	int mesh_node_count;
	int mesh_element_count;
	std::vector<Node> nodes;
	std::vector<Element> elements;
	std::vector<Material> materials;
	std::vector<BoundaryCondition> boundaries;
	std::vector<MaterialBlockRange> material_ranges;
	std::vector<FixedMaterialEntry> fixed_materials;

	std::vector<double> potentials;
	std::vector<double> current_sources;
	std::vector<bool> fixed_node;
	std::vector<double> fixed_value;

	// 疎行列: CSR
	std::vector<int> global_K_row_ptr;
	std::vector<int> global_K_col_idx;
	std::vector<double> global_K_values;

	// 誘電率行列 (E)
	std::vector<int> global_E_row_ptr;
	std::vector<int> global_E_col_idx;
	std::vector<double> global_E_values;

	// IC(0) preconditioner (lower triangular L in CSR and CSC for transpose)
	std::vector<int> L_row_ptr;
	std::vector<int> L_col_idx;
	std::vector<double> L_values;
	std::vector<int> L_csc_col_ptr;
	std::vector<int> L_csc_row_idx;
	std::vector<double> L_csc_values;
	std::vector<double> ic_scale;
	std::vector<double> preconditioner_diag;
	bool use_ssor_preconditioner = false;
	bool ic_built = false;
	// 組立途中の行リスト
	std::vector<std::vector<std::pair<int, double>>> global_K_rows;
	std::vector<std::vector<std::pair<int, double>>> global_E_rows;
	std::vector<double> global_F;

	double dt;
	int nstep;
	double unit;
	int final_step;
	double final_time;
	double final_max_diff;
	bool converged_early;
	int last_solver_iterations;
	double last_solver_time;
	double last_solver_relative_residual;
	bool last_solver_converged;
	std::string last_solver_method;
	std::vector<int> step_solver_iterations;
	std::vector<double> step_solver_times;
	std::string step_output_prefix;
	std::string report_output_path;
	std::string material_file_path;

	void add_fixed_node(int node, double value, bool overwrite);
	void build_sparse_from_rows(const std::vector<std::vector<std::pair<int, double>>> &rows,
		std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &values) const;
	void build_sparse_matrix();
	void build_time_step_matrix(double scale_eps,
		std::vector<int> &row_ptr, std::vector<int> &col_idx, std::vector<double> &values) const;
	void apply_dirichlet_constraints(std::vector<double> &rhs);
	void matvec(const std::vector<double> &x, std::vector<double> &y) const;
	void matvec_csr(const std::vector<int> &row_ptr, const std::vector<int> &col_idx,
		const std::vector<double> &values, const std::vector<double> &x, std::vector<double> &y) const;
	void write_analysis_report(const char *filename) const;
	// Incomplete Cholesky (IC(0)) preconditioner
	void build_incomplete_cholesky();
	void apply_ssor_preconditioner(const std::vector<double> &r, std::vector<double> &z) const;
	void apply_preconditioner(const std::vector<double> &r, std::vector<double> &z) const;

public:
	ElectrostaticAnalyzer();

	void set_output_paths(const std::string &step_prefix, const std::string &report_path);

	void read_mesh(const char *filename);                // in.dat
	void read_material_and_bc(const char *filename);     // sin.dat / thermo.dat
	void read_initial_potential(const char *filename);   // tmate.dat
	void read_analysis_config(const char *filename);     // sina.dat
	void write_potential_distribution(const char *filename); // tempa.dat001

	void solve();
	void apply_boundary_conditions();
	void assemble_global_matrix();
	void solve_linear_system();

	int get_num_nodes() const { return mesh_node_count; }
	int get_num_elements() const { return mesh_element_count; }
	int get_boundary_count() const { return (int)boundaries.size(); }
	int get_fixed_material_count() const { return (int)fixed_materials.size(); }
	int get_material_range_count() const { return (int)material_ranges.size(); }
	int get_material_count() const { return (int)materials.size(); }
	int get_final_step() const { return final_step; }
	double get_final_time() const { return final_time; }
	double get_final_max_diff() const { return final_max_diff; }
	bool get_converged_early() const { return converged_early; }
	int get_last_solver_iterations() const { return last_solver_iterations; }
	double get_last_solver_time() const { return last_solver_time; }
	const std::string &get_last_solver_method() const { return last_solver_method; }
	const std::vector<int> &get_step_solver_iterations() const { return step_solver_iterations; }
	const std::vector<double> &get_step_solver_times() const { return step_solver_times; }
	const std::vector<BoundaryCondition> &get_boundaries() const { return boundaries; }
	const std::vector<MaterialBlockRange> &get_material_ranges() const { return material_ranges; }
	const std::vector<FixedMaterialEntry> &get_fixed_materials() const { return fixed_materials; }
	const std::vector<Material> &get_materials() const { return materials; }
};

#endif
