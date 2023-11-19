
#include "R_stors.h"
#include "stors.h"

// return pdf value for the target dist
double f(double x, SEXP Rf, SEXP Renv)
{
    SEXP fcall, result, foo, val;
    double res;

    val = PROTECT(allocVector(REALSXP, 1));

    REAL(val)
    [0] = x;

    PROTECT(fcall = Rf_lang2(Rf, val));
    PROTECT(result = eval(fcall, Renv));
    PROTECT(foo = coerceVector(result, REALSXP));
    res = REAL(foo)[0];
    UNPROTECT(4);
    return (res);
}


SEXP print_cached_grids(void){
  
   Rprintf(" ========== \n");
   Rprintf(" x\n");
   for( int j=0; j < grids.incache; j++){
    for( int i = 0; i < (grids.grid[j].steps_number + 1) ; i++){
      Rprintf(" %f \n",grids.grid[j].x[i]);
    }
   }
   Rprintf(" ========== \n");
  
  R_RETURN_NULL
}
/*
SEXP pre_fetch( SEXP Rgrid , SEXP Rsize){
  
  double *obj = REAL(Rgrid);
  int size = asInteger(Rsize);

  for(int i =0; i < size; i++){
    
     __builtin_prefetch(obj + (sizeof(double) * i * N), 0, 3);
    
  }
  
  return (R_NilValue);
}
*/

SEXP stors(SEXP s_size, SEXP Rx, SEXP Rs_upper_lower, SEXP Rp_a, SEXP Rm,
           SEXP Rnormalized_areas, SEXP Runif_s, SEXP Rs_upper, SEXP Rlts, SEXP Rrts, SEXP Rf, SEXP Renv)
{
  int j, sample_size = asInteger(s_size), m = asInteger(Rm);
  
  double h_upper, u, u1, sample, f_sample, unif_s = asReal(Runif_s);
  
  double *upper = REAL(Rs_upper), *x = REAL(Rx), *s_upper_lower = REAL(Rs_upper_lower), *pa = REAL(Rp_a), *n_areas = REAL(Rnormalized_areas), *lts = REAL(Rlts), *rts = REAL(Rrts);
  
  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));
  
  double *results = REAL(Rresults);
  
  GetRNGstate();
  
  u1 = unif_rand();
  
  for (int i = 0; i < sample_size;)
  {
    
    if (u1 < n_areas[0])
    {
      
      sample = x[0] + (log(lts[0] + u1 * lts[1]) - lts[2]) * lts[3];
      
      h_upper = lts[4] * (sample - x[0]) + lts[2];
      
      u = unif_rand();
      
      if (u < f(sample, Rf, Renv) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }
      
      u1 = unif_rand();
    }
    else if (u1 > n_areas[1])
    {
      
      sample = x[m] + log1p((u1 * rts[0] - rts[1]) * rts[2]) * rts[3];
      
      h_upper = rts[4] * (sample - x[m]) + rts[5];
      
      u = unif_rand();
      
      if (u < f(sample, Rf, Renv) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }
      
      u1 = unif_rand();
    }
    else
    {
      
      u1 = (u1 - n_areas[0]) * unif_s;
      
      u1 *= m;
      
      j = (int)u1; // floor(u * m)
      
      u1 -= j;
      
      if (u1 < pa[j])
      {
        u1 = u1 * s_upper_lower[j];
        
        sample = x[j] + u1 * (x[j + 1] - x[j]);
        
        results[i] = sample;
        
        i++;
        
        if (i < sample_size)
        {
          u1 = unif_rand();
        }
      }
      else
      {
        
        double u0 = unif_rand();
        
        sample = x[j] + u0 * (x[j + 1] - x[j]);
        
        f_sample = f(sample, Rf, Renv);
        
        double uf = f_sample / upper[j];
        
        if (u1 < uf)
        {
          
          results[i] = sample;
          i++;
        }
        
        u1 = unif_rand();
      }
    }
  }
  
  PutRNGstate();
  
  UNPROTECT(1);
  
  return (Rresults);
}


