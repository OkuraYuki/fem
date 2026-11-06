// Build: g++ -O2 -std=c++17 tests/physics_regression.cpp src/electrostatic_analyzer.cpp -o /tmp/fem-physics-test
// Run with an existing temporary output directory as the only argument.
#include <vector>
#include <utility>
#include <string>
#define private public
#include "../src/electrostatic_analyzer.h"
#undef private
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>

static void check(double actual, double expected, const char *name) {
    if (!std::isfinite(actual) || std::fabs(actual - expected) > 1e-9) {
        std::fprintf(stderr, "FAIL %s: actual=%.15g expected=%.15g\n", name, actual, expected);
        std::exit(1);
    }
}

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    // Affine fields on a skew tetrahedron, including all local node orderings.
    int order[4] = {0,1,2,3};
    do {
        ElectrostaticAnalyzer a;
        a.nodes = {{0,0,0},{1,0,0},{1,1,0},{0,0,1}};
        a.elements = {Element{{order[0],order[1],order[2],order[3]},0}};
        a.materials = {Material{2,3,5,7,11,13,0,0}};
        a.assemble_global_matrix();
        for (int axis=0; axis<3; ++axis) {
            std::vector<double> field(4), product;
            for (int i=0; i<4; ++i)
                field[i] = axis==0 ? a.nodes[i].x : axis==1 ? a.nodes[i].y : a.nodes[i].z;
            a.matvec_csr(a.global_K_row_ptr,a.global_K_col_idx,a.global_K_values,field,product);
            double energy=0;
            for (int i=0; i<4; ++i) energy += field[i]*product[i];
            check(energy,(axis==0?2:axis==1?3:5)/6.,"conductive affine energy");
            a.matvec_csr(a.global_E_row_ptr,a.global_E_col_idx,a.global_E_values,field,product);
            energy=0;
            for (int i=0; i<4; ++i) energy += field[i]*product[i];
            check(energy,(axis==0?7:axis==1?11:13)/6.,"dielectric affine energy");
        }
    } while (std::next_permutation(order,order+4));

    // In a homogeneous proportional sigma/epsilon medium, charge-free potential
    // follows the imposed waveform exactly: the free node here is V_electrode/3.
    for (int mode=0; mode<7; ++mode) {
        ElectrostaticAnalyzer a;
        a.nodes={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        a.elements={Element{{0,1,2,3},0}};
        a.mesh_node_count=4; a.mesh_element_count=1;
        a.materials={Material{.01,.01,.01,1e-10,1e-10,1e-10,0,0}};
        a.fixed_node={false,true,true,true}; a.fixed_value={0,0,0,0};
        PulseVoltageEntry pulse{0,0,1,0,3e-9,0,100e-9,0,{1}};
        if(mode==1) pulse.phase=1e-9;
        if(mode==2) { pulse.t_rise=2e-9; pulse.t_fall=2e-9; }
        if(mode==3) { pulse.v_low=.2; pulse.v_high=1.2; }
        if(mode==4) { a.materials[0].sy=.02; a.materials[0].sz=.03; }
        a.pulse_voltages={pulse}; a.dt=1e-9; a.nstep=10;
        a.boundaries={{0,1,5,0}}; // Natural boundary must leave node 0 free.
        if(mode==5) {
            std::string path=std::string(argv[1])+"/time_config.dat";
            std::ofstream config(path);
            config << "# time parser regression\n0 0 0\n10 0 0 # nstep\n1e-9 .5e-9 0 # dt time0\n0 0\n";
            config.close();
            a.read_analysis_config(path.c_str());
            check(a.time0*1e9,.5,"parsed start time ns");
        }
        if(mode==6) { a.pulse_voltages.clear(); a.fixed_value[1]=1; }
        auto prefix=std::string(argv[1])+"/"+std::string(80,'x')+std::to_string(mode);
        a.set_output_paths(prefix+"_",prefix+".txt");
        a.voltage_difference_monitors={VoltageDifferenceMonitor{2,0,"free_node",{}, {}}};
        a.solve();
        if(mode==6) {
            check(a.voltage_difference_monitors[0].voltage_b[0],1./3.,"constant-voltage initialization");
            if(a.get_final_step()!=1 || !a.get_converged_early()) return 1;
            continue;
        }
        if(a.get_final_step()!=10 || a.get_converged_early()) return 1;
        if(!std::ifstream(prefix+"_001")) return 1;
        check(a.get_final_time()*1e9,(a.time0+10*a.dt)*1e9,"absolute final time ns");
        double previous_g=a.evaluate_pulse_voltage(pulse,a.time0);
        double expected=previous_g/3.;
        for(int step=1; step<=10; ++step) {
            double g=a.evaluate_pulse_voltage(pulse,a.time0+step*a.dt);
            if(mode==4)
                expected=(.3*expected+.01*g+.1*(g-previous_g))/.36;
            else expected=g/3.;
            check(a.voltage_difference_monitors[0].voltage_b[step-1],expected,"pulse response");
            previous_g=g;
        }
        if(mode==0) {
            const auto before=a.voltage_difference_monitors[0].voltage_b;
            a.solve();
            for(size_t i=0;i<before.size();++i)
                check(a.voltage_difference_monitors[0].voltage_b[i],before[i],"repeat solve node order");
        }
    }
    std::puts("PASS: affine fields (24 orderings, 3 axes, sigma and epsilon), 6 pulse cases plus constant-voltage initialization, natural boundary, nonzero start time, long output paths");
}
