#ifndef ELECTROSTATIC_ANALYZER_H
#define ELECTROSTATIC_ANALYZER_H

#include <vector>

struct Node {
	double x, y, z;
};

struct Element {
	int nodes[4];      // 4節点四面体
	int material;      // 0-based material id
};

struct Material {
	double sx, sy, sz; // 導電率
	double ro, cv;     // 互換用ダミー
};

struct BoundaryCondition {
	int node1, node2;  // 0-based node id
	int type;          // 1/4: Dirichlet, 5: Neumann(自然境界)
	double value;
};

class ElectrostaticAnalyzer {
private:
	std::vector<Node> nodes;
	std::vector<Element> elements;
	std::vector<Material> materials;
	std::vector<BoundaryCondition> boundaries;

	std::vector<double> potentials;
	std::vector<double> current_sources;
	std::vector<bool> fixed_node;
	std::vector<double> fixed_value;

	std::vector<std::vector<double>> global_K;
	std::vector<double> global_F;

	double dt;
	int nstep;
	double unit;

	void add_fixed_node(int node, double value, bool overwrite);

public:
	ElectrostaticAnalyzer();

	void read_mesh(const char *filename);                // in.dat
	void read_material_and_bc(const char *filename);     // sin.dat / thermo.dat
	void read_initial_potential(const char *filename);   // tmate.dat
	void read_analysis_config(const char *filename);     // sina.dat
	void write_potential_distribution(const char *filename); // tempa.dat001

	void solve();
	void apply_boundary_conditions();
	void assemble_global_matrix();
	void solve_linear_system();

	int get_num_nodes() const { return (int)nodes.size(); }
	int get_num_elements() const { return (int)elements.size(); }
};

#endif
