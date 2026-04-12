#include "pdesolver2d.h"

class SheetConductorSolver : public PDESolver2D
{
public:
  SheetConductorSolver(int n)
      : PDESolver2D(n)
      { }
  
  PDESolver2D *create_instance(int n) 
      { return new SheetConductorSolver(n); }
  
  double calc_new_value(double c, double t, double r, double b, double l, int p, int q);
};

double SheetConductorSolver::calc_new_value(double c, double t, double r, double b, double l, int p, int q)
{
  return 0.25 * (t + r + b + l);  // Laplace
}

////////////////////////////////////////////////////////

double sqr(double v) { return v*v; }

void set_boundary(PDESolver2D *solver)
{
  for(int p = 0; p < solver->size(); p++){
    for(int q = 0; q < solver->size(); q++){

      // outer boundary
      if( p==0 || q==0 || p==solver->size()-1 || q==solver->size()-1 )
        solver->pixel(p, q, NAN); // free edge

      // current injection
      double x = 2.0 * p / ( solver->size() - 1 ) - 1.0;
      double y = 2.0 * q / ( solver->size() - 1 ) - 1.0;

      // circular electrodes
      if ( sqr(x+0.5) + sqr(y) < sqr(0.05) ) {
        solver->mask(p, q, 1);    // fixed edge
        solver->pixel(p, q, -1.0);
      }
      if ( sqr(x-0.5) + sqr(y) < sqr(0.05) ) {
        solver->mask(p, q, 1);    // fixed edge
        solver->pixel(p, q, +1.0);
      }
    }
  }
}

PDESolver2D *solver = NULL;

// ^C が押されたら現在の計算を中断する
void signal_received(int sig)
{
  if(solver)
    solver->interrupt();
}

int main(int argc, const char *argv[])
{
  int initial_size = 17;
  
  solver = new SheetConductorSolver(initial_size);
  
  // ^C が押されたら signal_received へ飛ぶ
  signal(SIGINT, signal_received);
  
  // 徐々に解像度を上げながら計算する
  char fname[256];
  double delta = 1e-6;
  FILE *log;
  //  0   1   2    3    4    5     6
  // 17, 33, 65, 129, 257, 513, 1025
  for(int i = 0; i < 7; i++) {
    solver->set_speed(1);

    set_boundary(solver);
    sprintf(fname, "solver%dbefore.data", i);
    solver->dump(fname);

    sprintf(fname, "solver%d.log", i);
    log = fopen(fname, "w");
    int success = solver->solve(delta, log);
    fclose(log);

    if(success){
      sprintf(fname, "solver%dafter.data", i);
    } else {
      sprintf(fname, "solver%dinterrupted.data", i);
    }

    solver->dump(fname);

    fprintf(stderr, "\n*********** writing %s \n", fname);

    int delete_original = 1;
    solver = solver->create_double(delete_original);
  }
}