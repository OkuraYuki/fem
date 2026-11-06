#include <vector>
#include <utility>
#include <string>
#define private public
#include "../src/electrostatic_analyzer.h"
#undef private
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>

static void near(double a,double b,double tolerance,const char *name){
 if(!std::isfinite(a)||std::fabs(a-b)>tolerance){std::fprintf(stderr,"FAIL %s: %.17g expected %.17g\n",name,a,b);std::exit(1);}
}
static ElectrostaticAnalyzer tetra(const std::string &dir,int mode){
 ElectrostaticAnalyzer a;a.nodes={{0,0,0},{1,0,0},{0,1,0},{0,0,1}};a.elements={Element{{0,1,2,3},0}};
 a.mesh_node_count=4;a.mesh_element_count=1;a.materials={Material{.01,.01,.01,1e-10,1e-10,1e-10,0,0}};
 a.fixed_node={false,false,true,true};a.fixed_value={0,0,0,0};a.dt=1e-9;a.nstep=20;
 a.pulse_voltages={PulseVoltageEntry{1,0,1,0,(mode==1||mode==4)?100e-9:3e-9,mode==5?4e-9:0,100e-9,0,{1}}};
 a.circuit.enabled=true;a.circuit.material_id=1;a.circuit.nodes={1};a.circuit.coordinate_scale=mode==3?.001:1.;
 a.circuit.disconnect_time=(mode==1||mode==4||mode==5)?std::numeric_limits<double>::infinity():3e-9;
 a.circuit.series_resistance=(mode==1||mode==4||mode==5)?900:0;a.circuit.off_resistance=mode==2?std::numeric_limits<double>::infinity():900;
 a.circuit.external_capacitance=mode==0?1e-11:0;a.circuit.steady_initial_state=mode!=1;
 a.set_output_paths(dir+"/c"+std::to_string(mode)+"_",dir+"/c"+std::to_string(mode)+".txt");
 a.voltage_difference_monitors={{2,0,"interior",{},{}},{2,1,"electrode",{},{}}};return a;
}
int main(int argc,char**argv){
 if(argc!=2)return 2;
 std::string dir=argv[1];
 for(int mode=0;mode<6;++mode){auto a=tetra(dir,mode);a.solve();
  double G=.01/9*a.circuit.coordinate_scale,C=1e-10/9*a.circuit.coordinate_scale+a.circuit.external_capacitance;
  double expected=mode==1?0:(mode==4||mode==5)?.5:1;
  for(int step=0;step<=20;++step){
   if(step>0){
    if(mode==1||mode==4||mode==5) expected=(C/a.dt*expected+a.evaluate_pulse_voltage(a.pulse_voltages[0],step*a.dt)/900)/(C/a.dt+G+1./900);
    else if(step>3) expected=expected/(1+a.dt*(G+(mode==2?0:1./900))/C);
   }
   near(a.circuit.voltages[step],expected,2e-10,"RC electrode voltage");
   near(a.circuit.charges[step],C*expected,1e-18,"stored electrode charge");
   if(step>0 && ((mode==1||mode==4||mode==5)||step>3)) {
    double voltage=a.circuit.voltages[step];
    double current=(mode==1||mode==4||mode==5)?(a.circuit.source_voltages[step]-voltage)/900:(mode==2?0:-voltage/900);
    near((a.circuit.charges[step]-a.circuit.charges[step-1])/a.dt+G*voltage-current,0,1e-11,"electrode KCL");
   }
   if(step>0)near(a.voltage_difference_monitors[0].voltage_b[step-1],expected/3,2e-10,"interior field follows electrode");
  }
  // The capacitor voltage is continuous at the disconnection instant.
  if(mode!=1&&mode!=4&&mode!=5)near(a.circuit.voltages[3],a.circuit.voltages[2],2e-10,"switch continuity");
  if(mode==0){auto before=a.circuit.voltages;a.solve();for(size_t i=0;i<before.size();++i)near(a.circuit.voltages[i],before[i],2e-10,"repeat coupled solve");}
 }
 // Small extruded triangular prism, with z anisotropy and an imposed ramp.
 ElectrostaticAnalyzer a;a.nodes={{0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}};
 a.elements={Element{{0,1,2,3},0},Element{{1,2,3,4},0},Element{{2,3,4,5},0}};
 a.mesh_node_count=6;a.mesh_element_count=3;a.materials={Material{.01,.01,.008,1e-10,1e-10,1e-10,0,0}};
 a.fixed_node={false,true,true,false,true,true};a.fixed_value={0,1,0,0,1,0};
 a.pulse_voltages={PulseVoltageEntry{1,0,1,0,2e-9,4e-9,100e-9,0,{1,4}}};a.axial_symmetry=true;a.dt=1e-9;a.nstep=10;
 a.set_output_paths(dir+"/sym_",dir+"/sym.txt");a.voltage_difference_monitors={{2,0,"bottom",{},{}},{5,3,"top",{},{}}};a.solve();
 if(a.last_dof_count!=3||a.nodes.size()!=6||a.potentials.size()!=6)return 1;
 for(int step=1;step<=10;++step){double expected=.5*a.evaluate_pulse_voltage(a.pulse_voltages[0],step*a.dt);
  near(a.voltage_difference_monitors[0].voltage_b[step-1],expected,2e-10,"symmetric ramp");near(a.voltage_difference_monitors[1].voltage_b[step-1],expected,2e-10,"equal axial field");}
 auto before=a.voltage_difference_monitors[0].voltage_b;a.solve();for(size_t i=0;i<before.size();++i)near(a.voltage_difference_monitors[0].voltage_b[i],before[i],2e-10,"repeat symmetry solve");
 std::ifstream in(dir+"/sym_001");double v;int count=0;while(in>>v)++count;if(count!=6)return 1;
 std::puts("PASS: finite-R charging and steady initialization, finite-R falling source, resistive/open discharging, capacitance and SI scaling, switch continuity, charge and interior voltage, repeated solves, symmetric ramp and full-node output");
}
