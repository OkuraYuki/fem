#include <vector>
#include <utility>
#include <string>
#define private public
#include "../src/electrostatic_analyzer.h"
#undef private
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <numeric>

static void save(const std::string &p,const std::vector<double>&x){std::ofstream f(p,std::ios::binary);f.write(reinterpret_cast<const char*>(x.data()),x.size()*sizeof(double));}
static void collapse(ElectrostaticAnalyzer &a){
 const int n=2409,old=a.nodes.size();a.nodes.resize(n);
 auto matrix=[&](std::vector<int>&rp,std::vector<int>&ci,std::vector<double>&v){
  std::vector<std::vector<std::pair<int,double>>> rows(n);
  for(int r=0;r<old;++r)for(int k=rp[r];k<rp[r+1];++k)rows[r%n].push_back({ci[k]%n,v[k]});
  a.build_sparse_from_rows(rows,rp,ci,v);
 };
 matrix(a.global_K_row_ptr,a.global_K_col_idx,a.global_K_values);matrix(a.global_E_row_ptr,a.global_E_col_idx,a.global_E_values);
 std::vector<bool> fixed(n,false);std::vector<double> values(n,0),F(n,0);
 for(int i=0;i<old;++i){F[i%n]+=a.global_F[i];if(a.fixed_node[i]){if(fixed[i%n]&&values[i%n]!=a.fixed_value[i])std::exit(3);fixed[i%n]=true;values[i%n]=a.fixed_value[i];}}
 a.nodes.resize(n);a.potentials.resize(n);a.fixed_node=fixed;a.fixed_value=values;a.global_F=F;
 for(auto &p:a.pulse_voltages){for(int &i:p.nodes)i%=n;std::sort(p.nodes.begin(),p.nodes.end());p.nodes.erase(std::unique(p.nodes.begin(),p.nodes.end()),p.nodes.end());}
}
int main(int argc,char**argv){
 if(argc!=3)return 2;std::string mode=argv[1],out=argv[2],src="models/model_entyu_1";
 ElectrostaticAnalyzer a;a.read_mesh((src+"/in.dat").c_str());a.read_material_and_bc((src+"/sin.dat").c_str());a.read_initial_potential((src+"/tmate.dat").c_str());a.read_pulse_config((src+"/pulse.dat").c_str());a.dt=1e-9;a.nstep=26;
 if(mode=="iso")a.materials[1].sz=.01;
 if(mode=="fine"||mode=="coarse"){
  const int oldslabs=150,newslabs=mode=="fine"?300:75,nn=2409,ee=13968;
  auto oldnodes=a.nodes;auto oldelements=a.elements;auto oldfixed=a.fixed_node;auto oldvalues=a.fixed_value;
  for(int l=0;l<oldslabs;++l)for(int j=0;j<ee;++j){auto &x=oldelements[l*ee+j];auto &b=oldelements[j];if(x.material!=b.material)std::exit(4);for(int q=0;q<4;++q)if(x.nodes[q]!=b.nodes[q]+l*nn)std::exit(5);}
  a.nodes.resize((newslabs+1)*nn);a.fixed_node.assign(a.nodes.size(),false);a.fixed_value.assign(a.nodes.size(),0);
  for(int l=0;l<=newslabs;++l)for(int j=0;j<nn;++j){a.nodes[l*nn+j]=oldnodes[j];a.nodes[l*nn+j].z=1500.*l/newslabs;a.fixed_node[l*nn+j]=oldfixed[j];a.fixed_value[l*nn+j]=oldvalues[j];}
  a.elements.resize(newslabs*ee);for(int l=0;l<newslabs;++l)for(int j=0;j<ee;++j){auto t=oldelements[j];for(int &id:t.nodes)id+=l*nn;a.elements[l*ee+j]=t;}
  for(auto &p:a.pulse_voltages){std::vector<int> base;for(int i:p.nodes)if(i<nn)base.push_back(i);p.nodes.clear();for(int l=0;l<=newslabs;++l)for(int i:base)p.nodes.push_back(l*nn+i);}
 }
 a.potentials.assign(a.nodes.size(),0);a.apply_boundary_conditions();a.apply_pulse_voltages(0);a.assemble_global_matrix();
 if(mode=="collapsed")collapse(a);
 const auto kr=a.global_K_row_ptr,kc=a.global_K_col_idx;const auto kv=a.global_K_values,F=a.global_F;
 auto restore=[&](){a.global_K_row_ptr=kr;a.global_K_col_idx=kc;a.global_K_values=kv;a.global_F=F;};
 a.global_K_row_ptr=a.global_E_row_ptr;a.global_K_col_idx=a.global_E_col_idx;a.global_K_values=a.global_E_values;a.global_F.assign(a.nodes.size(),0);a.solve_linear_system();
 const auto initial=a.potentials;save(out+"/phi0.bin",initial);
 std::vector<double> ru,re;a.matvec_csr(kr,kc,kv,initial,ru);a.matvec_csr(a.global_E_row_ptr,a.global_E_col_idx,a.global_E_values,initial,re);
 double rmax=0,emax=0,zspan=0;for(size_t i=0;i<initial.size();++i){if(!a.fixed_node[i]){rmax=std::max(rmax,std::fabs(ru[i]));emax=std::max(emax,std::fabs(re[i]));}if(initial.size()>2409)zspan=std::max(zspan,std::fabs(initial[i]-initial[i%2409]));}
 std::ofstream summary(out+"/initial_metrics.txt");summary.precision(17);summary<<"nodes "<<initial.size()<<"\nKsigma_phi0_free_max "<<rmax<<"\nKepsilon_phi0_free_max "<<emax<<"\nmax_difference_from_bottom_layer "<<zspan<<"\n";
 // Lift the 2D reference back to the full mesh and inspect the free-node residual.
 if(mode!="collapsed"){
  std::ifstream in("output/poweroff_20261029/collapsed/phi0.bin",std::ios::binary);std::vector<double> ref(2409);in.read(reinterpret_cast<char*>(ref.data()),ref.size()*sizeof(double));if(!in)return 6;
  std::vector<double> lifted(initial.size());for(size_t i=0;i<lifted.size();++i)lifted[i]=ref[i%2409];
  a.matvec_csr(a.global_E_row_ptr,a.global_E_col_idx,a.global_E_values,lifted,re);
  double end=0,inside=0;int last=initial.size()/2409-1;
  for(size_t i=0;i<re.size();++i)if(!a.fixed_node[i]){double &r=(i/2409==0||int(i/2409)==last)?end:inside;r=std::max(r,std::fabs(re[i]));}
  summary<<"lifted_2D_E_residual_end_max "<<end<<"\nlifted_2D_E_residual_interior_max "<<inside<<"\n";
 }
 restore();a.potentials=initial;a.solve_linear_system();save(out+"/phi_steady.bin",a.potentials);double steadygap=0;for(size_t i=0;i<initial.size();++i)steadygap=std::max(steadygap,std::fabs(a.potentials[i]-initial[i]));summary<<"max_conductive_dielectric_field_difference "<<steadygap<<"\n";summary.close();
 a.potentials=initial;auto previous=initial;std::ofstream history(out+"/history.csv");history.precision(17);history<<"step,time_ns,min_V,max_V,bottom884_V,bottom828_V,center_r100_V\n";
 for(int step=1;step<=26;++step){restore();std::vector<int> rp,ci;std::vector<double> v;a.build_time_step_matrix(1./a.dt,rp,ci,v);a.global_K_row_ptr.swap(rp);a.global_K_col_idx.swap(ci);a.global_K_values.swap(v);std::vector<double> ep;a.matvec_csr(a.global_E_row_ptr,a.global_E_col_idx,a.global_E_values,previous,ep);for(size_t i=0;i<ep.size();++i)a.global_F[i]+=ep[i]/a.dt;a.potentials=previous;a.apply_pulse_voltages(step*a.dt);a.solve_linear_system();previous=a.potentials;
  auto range=std::minmax_element(previous.begin(),previous.end());int central=(int(previous.size()/2409)/2)*2409+824;
  history<<step<<","<<step<<","<<*range.first<<","<<*range.second<<","<<previous[883]<<","<<previous[827]<<","<<previous[central]<<"\n";history.flush();
  if(step==24||step==25)save(out+"/phi"+std::to_string(step)+".bin",previous);
 }
}
