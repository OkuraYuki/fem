#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cmath>

struct Node {
    double x, y, z;
};

struct Element {
    int nodes[4];      // 4節点四面体要素
    int material;
};

struct Material {
    double sx, sy, sz;  // 導電率 σ (thermal では κ)
    double ro, cv;      // ダミー（互換性維持）
};

struct BoundaryCondition {
    int node1, node2;
    int type;           // 0: Dirichlet (固定電位)
                        // 1: 全領域固定
                        // 4: 標準固定値
                        // 5: Neumann (自由端、電位勾配ゼロ)
    double value;       // 電位値 φ (Dirichlet用) または勾配値 (Neumann用)
};

class ElectrostaticAnalyzer {
private:
    std::vector<Node> nodes;
    std::vector<Element> elements;
    std::vector<Material> materials;
    std::vector<BoundaryCondition> boundaries;
    std::vector<double> potentials;     // 各節点の電位
    std::vector<bool> fixed_node;        // 固定電位ノードフラグ
    std::vector<double> fixed_value;     // 固定電位値
    std::vector<double> current_sources; // 電流源 I

    double dt;          // 時間刻幅（定常解析では不要だが互換性のため）
    int nstep;          // 反復回数
    double unit;        // 長さの単位

    // FEM用
    std::vector<std::vector<double>> global_K; // 大域剛性行列
    std::vector<double> global_F;              // 大域荷重ベクトル

    // 形状関数とその導関数
    double shape_function(int i, double xi, double eta, double zeta);
    void shape_derivatives(int i, double dN[3]);

    // Gauss積分点
    static const int ngauss = 4;
    double gauss_points[4][4]; // [point][xi,eta,zeta,weight]

public:
    ElectrostaticAnalyzer();
    // ファイル I/O
    void read_mesh(const char *filename);          // in.dat
    void read_material_and_bc(const char *filename); // thermo.dat
    void read_initial_potential(const char *filename); // tmate.dat
    void read_analysis_config(const char *filename); // sina.dat
    void write_potential_distribution(const char *filename); // tempa.dat001

    // 計算
    void solve();                      // Poisson方程式を解く
    void apply_boundary_conditions();
    void assemble_global_matrix();     // 大域剛性行列組立
    void solve_linear_system();        // 線形方程式を解く（反復法）

    // ユーティリティ
    int get_num_nodes() { return nodes.size(); }
    int get_num_elements() { return elements.size(); }
    void identify_boundary_nodes();  // 境界節点を特定
    std::vector<bool> is_boundary_node;  // 境界節点フラグ
};