#include "pdesolver2d.h"

PDESolver2D::PDESolver2D(int n)
{
  speed = 1;
  this->n = n;
  interrupted = 0;
  pixels = new double[n*n];
  masks = new int[n*n];
  clear();
}

PDESolver2D::~PDESolver2D()
{
  delete pixels;
  delete masks;
}

int PDESolver2D::solve(double delta, FILE *log)
{
  time_t started_at = time(NULL);
  
  int i;
  double diff = -1;
  int how_often = 53;
  for(i = 1; !interrupted; i++) {
    diff = step();
    if(log && (i % how_often == 0))
      log_step(i, diff, started_at, log);
    if(diff < delta) break;
  }
  if(log) log_step(i, diff, started_at, log);

  return !interrupted;
}

void PDESolver2D::log_step(int i, double diff, time_t started_at, FILE *log)
{
  if(log!=stderr)
    fprintf(log, "step (%6d, %e, %5.1fs)\n", i, diff, difftime(time(NULL), started_at));
  fprintf(stderr, "step (%6d, %e, %5.1fs)\r", i, diff, difftime(time(NULL), started_at));
}

double PDESolver2D::half_step(int even_or_odd)
{
  int p, q, nm1 = n-1;
  double result = 0;
  for(q = 1; q < nm1; q += 1 ){
    for(p = 1 + (even_or_odd  ^ (q % 2)); p < nm1; p += 2 ){
      if(mask(p, q)) continue;    // 固定端
      double value = pixel(p, q);
      if(isnan(value)) continue;  // 自由端
      double new_value = calc_new_value(
                  value,
                  replace_nan(p  , q-1, value),
                  replace_nan(p+1, q  , value),
                  replace_nan(p  , q+1, value),
                  replace_nan(p-1, q  , value),
                  p, q
              );
      pixel(p, q, value + speed * (new_value-value));
      
      double diff = fabs(new_value-value)/(fabs(new_value)+fabs(value));
      if( !isnan(diff) && diff > result ) result = diff;
    }
  }
  return result;
}

double PDESolver2D::step()
{
  return std::max( half_step(0), half_step(1) );
}

PDESolver2D *PDESolver2D::create_double(int delete_this)
{
  int p, q;
  PDESolver2D *result = create_instance(2*(n-1)+1);
  for(p = 0; p < result->size(); p++){
    for(q = 0; q < result->size(); q++){
      result->pixel(p, q, interpolate(p, q));
    }
  }
  if(delete_this) delete this;
  return result;
}

double PDESolver2D::interpolate(int p, int q)
{
  double i = 0;
  double result = 0, pix;
  int p2 = p/2, q2 = q/2;
  if(p % 2){  // need interpolation for p
    if(q % 2){  // need interpolation for q
      if (!isnan(pix = pixel(p2  , q2  )) ) { result += pix; i++; }
      if (!isnan(pix = pixel(p2+1, q2  )) ) { result += pix; i++; }
      if (!isnan(pix = pixel(p2  , q2+1)) ) { result += pix; i++; }
      if (!isnan(pix = pixel(p2+1, q2+1)) ) { result += pix; i++; }
      return i > 0 ? result / i : NAN;
    } else {    // not need interpolation for q
      if (!isnan(pix = pixel(p2  , q2  )) ) { result += pix; i++; }
      if (!isnan(pix = pixel(p2+1, q2  )) ) { result += pix; i++; }
      return i > 0 ? result / i : NAN;
    }
  } else {    // not need interpolation for p
    if(q % 2){  // need interpolation for q
      if (!isnan(pix = pixel(p2  , q2  )) ) { result += pix; i++; }
      if (!isnan(pix = pixel(p2  , q2+1)) ) { result += pix; i++; }
      return i > 0 ? result / i : NAN;
    } else {    // not need interpolation for q
      return pixel(p2, q2);
    }
  }
}

void PDESolver2D::dump(const char *file_name)
{
  FILE *f = fopen(file_name, "w");
  fwrite(pixels, sizeof(double), n*n, f);
  fclose(f);
}