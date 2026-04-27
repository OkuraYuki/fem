#include "electrostatic_analyzer.h"

#include <cstdio>
#include <cmath>
#include <algorithm>

ElectrostaticAnalyzer::ElectrostaticAnalyzer()
	: dt(0.0), nstep(200), unit(1.0)
{
}

void ElectrostaticAnalyzer::add_fixed_node(int node, double value, bool overwrite)
{
	if (node < 0 || node >= (int)fixed_node.size()) return;
	if (!fixed_node[node] || overwrite) {
		fixed_node[node] = true;
		fixed_value[node] = value;
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

	boundaries.assign(ntblk, BoundaryCondition{0, 0, 0, 0.0});
	for (int i = 0; i < ntblk; ++i) {
		int node1 = 0, node2 = 0, type = 0;
		double value = 0.0;
		if (fscanf(fp, "%d %d %d %lf", &node1, &node2, &type, &value) != 4) {
			printf("Invalid boundary record at %d in %s\n", i + 1, filename);
			fclose(fp);
			return;
		}
		boundaries[i].node1 = node1 - 1; // 1-based -> 0-based
		boundaries[i].node2 = node2 - 1;
		boundaries[i].type = type;
		boundaries[i].value = value;
	}

	// 要素の材料番号割り当て
	for (int i = 0; i < ntmblk; ++i) {
		int start = 0, end = 0, mat = 0;
		if (fscanf(fp, "%d %d %d", &start, &end, &mat) != 3) {
			printf("Invalid material block record at %d in %s\n", i + 1, filename);
			fclose(fp);
			return;
		}
		int mat_id = mat - 1; // 1-based -> 0-based
		for (int e = start - 1; e <= end - 1; ++e) {
			if (0 <= e && e < (int)elements.size()) {
				elements[e].material = mat_id;
			}
		}
	}

	materials.assign(ntmate, Material{1.0, 1.0, 1.0, 0.0, 0.0});
	for (int i = 0; i < ntmate; ++i) {
		if (fscanf(fp, "%lf %lf %lf %lf", &materials[i].sx, &materials[i].sy, &materials[i].sz, &materials[i].ro) != 4) {
			// 5列形式にも対応
			fseek(fp, 0, SEEK_CUR);
			double cv = 0.0;
			if (fscanf(fp, "%lf %lf %lf %lf %lf", &materials[i].sx, &materials[i].sy, &materials[i].sz, &materials[i].ro, &cv) == 5) {
				materials[i].cv = cv;
			} else {
				printf("Invalid material record at %d in %s\n", i + 1, filename);
				fclose(fp);
				return;
			}
		}
	}

	int nqt = 0;
	if (fscanf(fp, "%d", &nqt) != 1) {
		nqt = 0;
	}

	current_sources.assign(nodes.size(), 0.0);
	for (int i = 0; i < nqt; ++i) {
		int elem_id = 0;
		double ix = 0.0, iy = 0.0, iz = 0.0;
		if (fscanf(fp, "%d %lf %lf %lf", &elem_id, &ix, &iy, &iz) != 4) break;
		int e = elem_id - 1;
		if (0 <= e && e < (int)elements.size()) {
			double per_node = (ix + iy + iz) / 4.0;
			for (int k = 0; k < 4; ++k) {
				int n = elements[e].nodes[k];
				if (0 <= n && n < (int)current_sources.size()) current_sources[n] += per_node;
			}
		}
	}

	fclose(fp);
	printf("Read material & BC: %d materials, %d boundary conditions, %d current sources\n", ntmate, ntblk, nqt);
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

	int flags[7] = {0};
	int params[7] = {0};
	double time0 = 0.0;
	double freq = 0.0;

	fscanf(fp, "%d %d %d %d %d %d %d", &flags[0], &flags[1], &flags[2], &flags[3], &flags[4], &flags[5], &flags[6]);
	fscanf(fp, "%d %d %d %d %d %d %d", &params[0], &params[1], &params[2], &params[3], &params[4], &params[5], &params[6]);
	fscanf(fp, "%lf %lf", &dt, &time0);
	fscanf(fp, "%lf", &freq);

	nstep = params[0] > 0 ? params[0] : 200;

	fclose(fp);
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
	global_K.assign(n, std::vector<double>(n, 0.0));
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
		double sigma = (materials[mat].sx + materials[mat].sy + materials[mat].sz) / 3.0;

		double volume = std::fabs(detJ) / 6.0;
		int gid[4] = {i0, i1, i2, i3};

		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b) {
				double dot = gradN[a][0] * gradN[b][0] + gradN[a][1] * gradN[b][1] + gradN[a][2] * gradN[b][2];
				global_K[gid[a]][gid[b]] += sigma * dot * volume;
			}
			global_F[gid[a]] += (gid[a] < (int)current_sources.size() ? current_sources[gid[a]] * volume / 4.0 : 0.0);
		}
	}

	printf("Assembled global matrix: %dx%d\n", n, n);
}

void ElectrostaticAnalyzer::solve_linear_system()
{
	const int n = (int)nodes.size();
	const int max_iter = std::max(10, nstep);
	const double tol = 1e-8;

	for (int iter = 0; iter < max_iter; ++iter) {
		double max_diff = 0.0;

		for (int i = 0; i < n; ++i) {
			if (fixed_node[i]) {
				potentials[i] = fixed_value[i];
				continue;
			}

			double aii = global_K[i][i];
			if (std::fabs(aii) < 1e-20) continue;

			double rhs = global_F[i];
			for (int j = 0; j < n; ++j) {
				if (j == i) continue;
				rhs -= global_K[i][j] * potentials[j];
			}

			double new_phi = rhs / aii;
			double diff = std::fabs(new_phi - potentials[i]);
			if (diff > max_diff) max_diff = diff;
			potentials[i] = new_phi;
		}

		if (max_diff < tol) {
			printf("Converged at iteration %d\n", iter + 1);
			break;
		}
	}

	printf("Solved linear system with Gauss-Seidel\n");
}

void ElectrostaticAnalyzer::solve()
{
	printf("Starting electrostatic analysis...\n");

	potentials.assign(nodes.size(), 0.0);
	if (fixed_node.size() != nodes.size()) {
		fixed_node.assign(nodes.size(), false);
		fixed_value.assign(nodes.size(), 0.0);
	}

	apply_boundary_conditions();
	assemble_global_matrix();
	solve_linear_system();

	printf("Analysis completed.\n");
}

void ElectrostaticAnalyzer::write_potential_distribution(const char *filename)
{
	FILE *fp = fopen(filename, "w");
	if (!fp) {
		printf("Cannot open %s\n", filename);
		return;
	}

	int col = 0;
	for (size_t i = 0; i < potentials.size(); ++i) {
		fprintf(fp, "%12.4E", potentials[i]);
		col++;
		if (col == 6 || i + 1 == potentials.size()) {
			fprintf(fp, "\n");
			col = 0;
		}
	}

	fclose(fp);
	printf("Wrote %s: %d nodal potentials\n", filename, (int)potentials.size());
}
