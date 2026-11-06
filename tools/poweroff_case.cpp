#include <vector>
#include <utility>
#include <string>
#define private public
#include "../src/electrostatic_analyzer.h"
#undef private
#include <cstdio>
#include <cstdlib>
// model_dir output_dir dt_ns end_ns rise_ns high_ns fall_ns phase_ns
int main(int argc,char **argv){
 if(argc!=9)return 2;
 ElectrostaticAnalyzer a;std::string src=argv[1],out=argv[2];
 a.read_mesh((src+"/in.dat").c_str());
 a.read_material_and_bc((src+"/sin.dat").c_str());
 a.read_initial_potential((src+"/tmate.dat").c_str());
 a.read_pulse_config((src+"/pulse.dat").c_str());
 a.read_analysis_config((src+"/sina.dat.t").c_str());
 a.dt=std::atof(argv[3])*1e-9;a.nstep=int(std::atof(argv[4])*1e-9/a.dt+.5);
 for(auto &p:a.pulse_voltages){p.t_rise=std::atof(argv[5])*1e-9;p.t_high=std::atof(argv[6])*1e-9;p.t_fall=std::atof(argv[7])*1e-9;p.t_low=1000e-9;p.phase=std::atof(argv[8])*1e-9;}
 for(int node:{181500,181254,180954,180733,180684})a.voltage_difference_monitors.push_back({180675,node-1,"r_node_"+std::to_string(node),{}, {}});
 a.set_output_paths(out+"/tempa.dat",out+"/analysis_summary.txt");a.solve();
 FILE *fp=fopen((out+"/timeseries.csv").c_str(),"w");
 fprintf(fp,"time_ns,input_V,r100_V,r500_V,r1200_V,r2000_V,r2400_V\n");
 for(int i=0;i<a.final_step;++i){double t=(i+1)*a.dt;fprintf(fp,"%.12g,%.12g",t*1e9,a.evaluate_pulse_voltage(a.pulse_voltages[0],t));for(auto &m:a.voltage_difference_monitors)fprintf(fp,",%.12g",m.voltage_b[i]);fprintf(fp,"\n");}fclose(fp);
}
