#include <vector>
#include <utility>
#include <string>
#define private public
#include "../src/electrostatic_analyzer.h"
#undef private
#include <fstream>
int main(){
 std::ifstream in("output/poweroff_20261029/fortran_kernel/samples.dat");std::ofstream out("output/poweroff_20261029/fortran_kernel/cpp_matrices.dat");out.precision(17);int samples;in>>samples;
 for(int n=0;n<samples;++n){ElectrostaticAnalyzer a;a.nodes.resize(4);for(auto &p:a.nodes)in>>p.x>>p.y>>p.z;a.elements={Element{{0,1,2,3},0}};a.materials={Material{.01,.01,.008,1e-10,1e-10,1e-10,0,0}};a.assemble_global_matrix();
 for(int r=0;r<4;++r)for(int c=0;c<4;++c){double k=0,e=0;for(int p=a.global_K_row_ptr[r];p<a.global_K_row_ptr[r+1];++p)if(a.global_K_col_idx[p]==c)k=a.global_K_values[p];for(int p=a.global_E_row_ptr[r];p<a.global_E_row_ptr[r+1];++p)if(a.global_E_col_idx[p]==c)e=a.global_E_values[p];out<<k<<" "<<e/1e-9<<"\n";}
 }
}
