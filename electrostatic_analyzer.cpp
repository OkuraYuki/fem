#include "electrostatic_analyzer.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

ElectrostaticAnalyzer::ElectrostaticAnalyzer() {
    // Gauss積分点の初期化 (4点積分)
    // 参考: TetrahedronのGauss積分点
    double a = (5.0 - sqrt(5.0)) / 20.0;
    double b = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
    double w = 1.0 / 24.0;

    gauss_points[0][0] = a; gauss_points[0][1] = a; gauss_points[0][2] = a; gauss_points[0][3] = w;
    gauss_points[1][0] = b; gauss_points[1][1] = a; gauss_points[1][2] = a; gauss_points[1][3] = w;
    gauss_points[2][0] = a; gauss_points[2][1] = b; gauss_points[2][2] = a; gauss_points[2][3] = w;
    gauss_points[3][0] = a; gauss_points[3][1] = a; gauss_points[3][2] = b; gauss_points[3][3] = w;
}

// 形状関数
double ElectrostaticAnalyzer::shape_function(int i, double xi, double eta, double zeta) {
    switch (i) {
        case 0: return 1.0 - xi - eta - zeta;
        case 1: return xi;
        case 2: return eta;
        case 3: return zeta;
        default: return 0.0;
    }
}

// 形状関数の導関数
void ElectrostaticAnalyzer::shape_derivatives(int i, double dN[3]) {
    switch (i) {
        case 0: dN[0] = -1.0; dN[1] = -1.0; dN[2] = -1.0; break;
        case 1: dN[0] =  1.0; dN[1] =  0.0; dN[2] =  0.0; break;
        case 2: dN[0] =  0.0; dN[1] =  1.0; dN[2] =  0.0; break;
        case 3: dN[0] =  0.0; dN[1] =  0.0; dN[2] =  1.0; break;
    }
}

// in.dat 読み込み（共通）
void ElectrostaticAnalyzer::read_mesh(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return;
    }

    int num_nodes, num_elements;
    int dummy1, dummy2;

    fscanf(fp, "%d %d %d %d %lf", &num_nodes, &num_elements,
           &dummy1, &dummy2, &unit);

    // 要素情報（1-based → 0-based）
    elements.resize(num_elements);
    for (int i = 0; i < num_elements; i++) {
        int n0, n1, n2, n3;
        fscanf(fp, "%d %d %d %d", &n0, &n1, &n2, &n3);
        elements[i].nodes[0] = n0 - 1;
        elements[i].nodes[1] = n1 - 1;
        elements[i].nodes[2] = n2 - 1;
        elements[i].nodes[3] = n3 - 1;
        elements[i].material = 0;
    }

    // 節点座標
    nodes.resize(num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        fscanf(fp, "%le %le %le", &nodes[i].x, &nodes[i].y, &nodes[i].z);
    }

    fclose(fp);

    printf("Read mesh: %d nodes, %d elements\n", num_nodes, num_elements);
}

// thermo.dat 読み込み（電位分布用に改良）
void ElectrostaticAnalyzer::read_material_and_bc(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return;
    }

    int ntblk, ntmblk, ntmate, npt2, nlv;

    fscanf(fp, "%d %d %d %d %d", &ntblk, &ntmblk, &ntmate, &npt2, &nlv);

    // 境界条件（1-based → 0-based、Dirichlet固定）
    boundaries.resize(ntblk);
    for (int i = 0; i < ntblk; i++) {
        int n1, n2;
        fscanf(fp, "%d %d %d %lf", &n1, &n2,
               &boundaries[i].type,
               &boundaries[i].value);  // 電位値 φ [V]
        boundaries[i].node1 = n1 - 1;
        boundaries[i].node2 = n2 - 1;
    }

    // 要素の材料番号（element id も 1-based → 0-based）
    for (int i = 0; i < ntmblk; i++) {
        int start, end, mat;
        fscanf(fp, "%d %d %d", &start, &end, &mat);
        start -= 1;
        end -= 1;
        int mat_id = mat - 1;  // 1-based material id → 0-based
        for (int j = start; j <= end; j++) {
            if (j >= 0 && j < elements.size()) {
                elements[j].material = mat_id;
            }
        }
    }

    // 材料データ（導電率）
    materials.resize(ntmate);
    for (int i = 0; i < ntmate; i++) {
        fscanf(fp, "%le %le %le %le %le",
               &materials[i].sx,      // σx
               &materials[i].sy,      // σy
               &materials[i].sz,      // σz
               &materials[i].ro,      // ダミー
               &materials[i].cv);     // ダミー
    }

    // 電流源データ
    int nqt;
    fscanf(fp, "%d", &nqt);
    current_sources.resize(nodes.size(), 0.0);

    for (int i = 0; i < nqt; i++) {
        int elem_id;
        double ix, iy, iz;
        fscanf(fp, "%d %le %le %le", &elem_id, &ix, &iy, &iz);
        int elem_index = elem_id - 1;
        // 電流源をノードに分配
        if (elem_index >= 0 && elem_index < elements.size()) {
            double current_per_node = (ix + iy + iz) / 4.0;
            for (int j = 0; j < 4; j++) {
                int node_id = elements[elem_index].nodes[j];
                if (node_id >= 0 && node_id < current_sources.size()) {
                    current_sources[node_id] += current_per_node;
                }
            }
        }
    }

    fclose(fp);

    printf("Read material & BC: %d materials, %d boundary conditions, "
           "%d current sources\n", ntmate, ntblk, nqt);
}

// sina.dat 読み込み（共通）
void ElectrostaticAnalyzer::read_analysis_config(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return;
    }

    int flags[7], params[7];
    fscanf(fp, "%d %d %d %d %d %d %d",
           &flags[0], &flags[1], &flags[2], &flags[3],
           &flags[4], &flags[5], &flags[6]);

    fscanf(fp, "%d %d %d %d %d %d %d",
           &params[0], &params[1], &params[2], &params[3],
           &params[4], &params[5], &params[6]);

    double time0;
    fscanf(fp, "%le %le", &dt, &time0);

    double freq;
    fscanf(fp, "%le", &freq);

    this->nstep = params[0];  // 反復回数

    fclose(fp);

    printf("Read config: dt=%.3e, nstep=%d\n", dt, nstep);
}

// 境界節点を特定
void ElectrostaticAnalyzer::identify_boundary_nodes()
{
    is_boundary_node.assign(nodes.size(), false);
    std::vector<int> node_element_count(nodes.size(), 0);

    // 各節点が属する要素数をカウント
    for (const auto &elem : elements) {
        for (int j = 0; j < 4; j++) {
            int node_id = elem.nodes[j];
            if (node_id < nodes.size()) {
                node_element_count[node_id]++;
            }
        }
    }

    // 要素数が1以下の節点を境界節点とする
    int boundary_count = 0;
    for (int i = 0; i < nodes.size(); i++) {
        if (node_element_count[i] <= 1) {
            is_boundary_node[i] = true;
            boundary_count++;
        }
    }

    printf("Identified %d boundary nodes\n", boundary_count);
}

void ElectrostaticAnalyzer::apply_boundary_conditions()
{
    for (const auto &bc : boundaries) {
        if (bc.type == 1 || bc.type == 4) {
            // 型1または4は両端ノードを固定 (Dirichlet)
            fixed_node[bc.node1] = true;
            fixed_value[bc.node1] = bc.value;
            fixed_node[bc.node2] = true;
            fixed_value[bc.node2] = bc.value;
        } else if (bc.type == 5) {
            // 型5: Neumann条件 (自由端、電位勾配ゼロ)
            // 境界節点を固定せず、自然境界条件を適用
            // 何もしない (デフォルトでNeumann)
        } else {
            // その他は片側ノードを固定 (Dirichlet)
            fixed_node[bc.node1] = true;
            fixed_value[bc.node1] = bc.value;
        }
    }

    // 固定電位を potentials に反映
    for (int i = 0; i < nodes.size(); ++i) {
        if (i < fixed_node.size() && fixed_node[i]) {
            potentials[i] = fixed_value[i];
        }
    }
}

void ElectrostaticAnalyzer::assemble_global_matrix()
{
    int n = nodes.size();
    global_K.assign(n, std::vector<double>(n, 0.0));
    global_F.assign(n, 0.0);

    for (size_t e = 0; e < elements.size(); ++e) {
        const Element &elem = elements[e];
        int mat_id = elem.material;
        if (mat_id >= materials.size()) continue;
        Material mat = materials[mat_id];

        // 要素ノードの座標
        Node n0 = nodes[elem.nodes[0]];
        Node n1 = nodes[elem.nodes[1]];
        Node n2 = nodes[elem.nodes[2]];
        Node n3 = nodes[elem.nodes[3]];

        // Jacobi行列の計算 (体積)
        double J[3][3] = {
            {n1.x - n0.x, n2.x - n0.x, n3.x - n0.x},
            {n1.y - n0.y, n2.y - n0.y, n3.y - n0.y},
            {n1.z - n0.z, n2.z - n0.z, n3.z - n0.z}
        };
        double detJ = J[0][0]*(J[1][1]*J[2][2] - J[1][2]*J[2][1]) -
                      J[0][1]*(J[1][0]*J[2][2] - J[1][2]*J[2][0]) +
                      J[0][2]*(J[1][0]*J[2][1] - J[1][1]*J[2][0]);
        if (detJ <= 0) continue; // 無効要素
        double volume = detJ / 6.0;

        // Jの逆行列
        double invJ[3][3];
        invJ[0][0] = (J[1][1]*J[2][2] - J[1][2]*J[2][1]) / detJ;
        invJ[0][1] = (J[0][2]*J[2][1] - J[0][1]*J[2][2]) / detJ;
        invJ[0][2] = (J[0][1]*J[1][2] - J[0][2]*J[1][1]) / detJ;
        invJ[1][0] = (J[1][2]*J[2][0] - J[1][0]*J[2][2]) / detJ;
        invJ[1][1] = (J[0][0]*J[2][2] - J[0][2]*J[2][0]) / detJ;
        invJ[1][2] = (J[0][2]*J[1][0] - J[0][0]*J[1][2]) / detJ;
        invJ[2][0] = (J[1][0]*J[2][1] - J[1][1]*J[2][0]) / detJ;
        invJ[2][1] = (J[0][1]*J[2][0] - J[0][0]*J[2][1]) / detJ;
        invJ[2][2] = (J[0][0]*J[1][1] - J[0][1]*J[1][0]) / detJ;

        // 局所剛性行列
        double Ke[4][4] = {0};

        for (int gp = 0; gp < ngauss; ++gp) {
            double xi = gauss_points[gp][0];
            double eta = gauss_points[gp][1];
            double zeta = gauss_points[gp][2];
            double weight = gauss_points[gp][3];

            // 物理座標での導関数 B = invJ * dN/dξ
            double B[4][3];
            for (int i = 0; i < 4; ++i) {
                double dN[3];
                shape_derivatives(i, dN);
                for (int j = 0; j < 3; ++j) {
                    B[i][j] = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        B[i][j] += invJ[j][k] * dN[k];
                    }
                }
            }

            // Ke += B^T * D * B * detJ * weight
            double sigma = mat.sx;
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    double sum = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        sum += B[i][k] * B[j][k];
                    }
                    Ke[i][j] += sigma * sum * detJ * weight;
                }
            }
        }

        // 大域行列に組み立て
        for (int i = 0; i < 4; ++i) {
            int gi = elem.nodes[i];
            for (int j = 0; j < 4; ++j) {
                int gj = elem.nodes[j];
                global_K[gi][gj] += Ke[i][j];
            }
            // 荷重ベクトル (電流源)
            global_F[gi] += current_sources[gi] * volume / 4.0; // 簡易分配
        }
    }

    printf("Assembled global matrix: %dx%d\n", n, n);
}

void ElectrostaticAnalyzer::solve_linear_system()
{
    int n = nodes.size();
    int max_iter = nstep; // sina.datから
    double tol = 1e-6;

    // Gauss-Seidel法
    for (int iter = 0; iter < max_iter; ++iter) {
        double max_diff = 0.0;
        for (int i = 0; i < n; ++i) {
            // 固定ポテンシャルを持つノードはスキップ
            if (i < fixed_node.size() && fixed_node[i]) continue;

            if (global_K[i][i] == 0.0) continue; // 対角要素0はスキップ

            double sum = global_F[i];
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum -= global_K[i][j] * potentials[j];
                }
            }
            double new_phi = sum / global_K[i][i];
            double diff = fabs(new_phi - potentials[i]);
            if (diff > max_diff) max_diff = diff;
            potentials[i] = new_phi;
        }
        if (max_diff < tol) {
            printf("Converged at iteration %d\n", iter);
            break;
        }
    }

    printf("Solved linear system with Gauss-Seidel\n");
}

// tmate.dat 読み込み（初期電位）
void ElectrostaticAnalyzer::read_initial_potential(const char *filename)
{
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return;
    }

    fixed_node.assign(nodes.size(), false);
    fixed_value.assign(nodes.size(), 0.0);

    int num_initial;
    fscanf(fp, "%d", &num_initial);

    for (int i = 0; i < num_initial; i++) {
        int mat_id;
        double value;
        fscanf(fp, "%d %le", &mat_id, &value);
        int material_index = mat_id - 1;

        // 材質番号 material_index の要素のノードに初期電位を設定
        for (const auto &elem : elements) {
            if (elem.material == material_index) {
                for (int j = 0; j < 4; j++) {
                    int node_id = elem.nodes[j];
                    if (node_id >= 0 && node_id < fixed_node.size()) {
                        fixed_node[node_id] = true;
                        fixed_value[node_id] = value;
                    }
                }
            }
        }
    }

    fclose(fp);

    printf("Read initial potential: %d entries\n", num_initial);
}

// tempa.dat001 出力
void ElectrostaticAnalyzer::write_potential_distribution(const char *filename)
{
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        fprintf(stderr, "Cannot open %s\n", filename);
        return;
    }

    int count = 0;
    for (int i = 0; i < potentials.size(); i++) {
        fprintf(fp, "%12.4E", potentials[i]);  // 電位値 φ [V]
        count++;

        if (count % 6 == 0 || i == potentials.size() - 1) {
            fprintf(fp, "\n");
            count = 0;
        }
    }

    fclose(fp);

    printf("Wrote %s: %d nodal potentials\n",
           filename, (int)potentials.size());
}

// 電位分布の計算（Poisson方程式）
void ElectrostaticAnalyzer::solve()
{
    printf("Starting electrostatic analysis...\n");

    int n = nodes.size();
    potentials.assign(n, 0.0);
    if (fixed_node.size() != n) {
        fixed_node.assign(n, false);
        fixed_value.assign(n, 0.0);
    }

    // 境界節点を特定
    identify_boundary_nodes();

    // 境界条件を適用
    apply_boundary_conditions();

    // 大域剛性行列を組立
    assemble_global_matrix();

    // 線形方程式を解く
    solve_linear_system();

    printf("Analysis completed.\n");
}