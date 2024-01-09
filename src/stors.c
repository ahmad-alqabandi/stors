
#include "stors.h"

#include "cache.h"


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

SEXP stors(SEXP s_size, SEXP R_Cnum, SEXP Rf, SEXP Renv)
{
  
  int j, sample_size = asInteger(s_size), Cnum = asInteger(R_Cnum);
  
  struct grid g = grids.grid[Cnum];
  
  if(g.exist == 0){
    REprintf("you need to optimize your destribution grid first");
    R_RETURN_NULL
  }
  
  double h_upper, u, u1, sample, f_sample;
  
  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));
  
  double *results = REAL(Rresults);
  
  GetRNGstate();
  
  u1 = unif_rand();
  
  for (int i = 0; i < sample_size;)
  {
    
    if (u1 < g.sampling_probabilities[0])
    {
      
      sample = g.x[0] + (log( g.lt_properties[0] + u1 * g.lt_properties[1]) - g.lt_properties[2]) * g.lt_properties[3];
      
      h_upper = g.lt_properties[4] * (sample - g.x[0]) + g.lt_properties[2];
      
      u = unif_rand();
      
      if (u < f(sample, Rf, Renv) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }
      
      u1 = unif_rand();
    }
    else if (u1 > g.sampling_probabilities[1])
    {
      
      sample = g.x[g.steps_number] + log1p((u1 * g.rt_properties[0] - g.rt_properties[1]) * g.rt_properties[2]) * g.rt_properties[3];
      
      h_upper = g.rt_properties[4] * (sample - g.x[g.steps_number]) + g.rt_properties[5];
      
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
      
      u1 = (u1 - g.sampling_probabilities[0]) * g.unif_scaler;
      
      u1 *= g.steps_number;
      
      j = (int)u1; // floor(u * m)
      
      u1 -= j;
      
      
      if (u1 < g.p_a[j])
      {
        u1 = u1 * g.s_upper_lower[j];
        
        sample = g.x[j] + u1 * (g.x[j + 1] - g.x[j]);
        
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
        
        sample = g.x[j] + u0 * (g.x[j + 1] - g.x[j]);
        
        f_sample = f(sample, Rf, Renv);
        
        double uf = f_sample /g.s_upper[j];
        
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



SEXP stors_trunc_nav(SEXP Rcnum ,SEXP Rlx, SEXP Rrx){
  
  int Cnum = asInteger(Rcnum);
  
  
  struct grid g = grids.grid[Cnum];
  
  if(g.exist == 0){
    REprintf("you need to optimize your destribution grid first");
    R_RETURN_NULL
  }
  
  int i=0;
  
  double xlr[2], total_area = g.areas[0]+ g.areas[1]+ g.areas[2], ixlr[2];
  

  double hu_x;
  

  
  xlr[0] = asReal(Rlx);
  xlr[1] = asReal(Rrx);
  
  SEXP Rresults = PROTECT((allocVector(REALSXP, 4)));
  
  double *results = REAL(Rresults);
  
  for( int j = 0; j < 2; j++){ 
    
    
    

    if(xlr[j] > g.x[g.steps_number]){
      
      hu_x =  g.rt_properties[4] * (xlr[j] - g.x[g.steps_number]) + g.rt_properties[5];
      
      results[j] =( g.areas[0] + g.areas[1] + (g.rt_properties[3] * (exp(hu_x) - exp(g.rt_properties[5]))) )/ total_area;
      
      ixlr[j] = -1;
      

    }else{
      
  while( i < g.steps_number + 1){
    
    if(xlr[j] < g.x[i]){
      
      

      if(i == 0){
        
        hu_x =  g.lt_properties[4] * (xlr[j] - g.x[0]) + g.lt_properties[2];
        
        results[j] = (g.lt_properties[3] * (exp(hu_x) - g.lt_properties[0])) / total_area;
        ixlr[j] = -1;

      }else{
        
  results[j] =(g.areas[0] + g.alpha * (i - 1) + (xlr[j] - g.x[i-1]) * g.s_upper[i-1]) / total_area;
  if( j == 0){
    ixlr[j] = i-1;
    
  }else{
    ixlr[j] = i;
    
  }
}
break;
    }
    
    i++;
    
  }
}


  }
  
  results[2] = ixlr[0];
  results[3] = ixlr[1];
  
  
  UNPROTECT(1);
  
  return (Rresults);
  
}





SEXP stors_trunc(SEXP s_size, SEXP Rcnum, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir,SEXP Rf, SEXP Renv){
  
  
  int j, sample_size = asInteger(s_size);
  double csl = asReal(Rcsl), csr = asReal(Rcsr), xl = asReal(Rxl), xr = asReal(Rxr), tmpxl, tempxr ;
  int il = asInteger(Ril), ir = asInteger(Rir);
  
  int Cnum = asInteger(Rcnum);
  
  
  struct grid g = grids.grid[Cnum];
  
  if(g.exist == 0){
    REprintf("you need to optimize your destribution grid first");
    R_RETURN_NULL
  }
  
  
  if(il != -1){
    tmpxl = g.x[il];
    g.x[il] = xl;
  }

  if(ir != -1){
    tempxr =g.x[ir];
    g.x[ir] = xr;
  }


  double h_upper, u;
  

  double  u1, sample, f_sample;
  
  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));
  
  double *results = REAL(Rresults);
  
  GetRNGstate();
  
  u1 = unif_rand();
  u1 = csl + u1 * (csr-csl);
  
  for (int i = 0; i < sample_size;)
  {
    
    

    if (u1 < g.sampling_probabilities[0])
    {

      sample = g.x[0] + (log( g.lt_properties[0] + u1 * g.lt_properties[1]) - g.lt_properties[2]) * g.lt_properties[3];
      h_upper = g.lt_properties[4] * (sample - g.x[0]) + g.lt_properties[2];
      u = unif_rand();
      if (u < f(sample, Rf, Renv) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }
      u1 = unif_rand();
      
      u1 = csl + u1 * (csr-csl);
      

      
    }else if(u1 > g.sampling_probabilities[1]){

        

        sample = g.x[g.steps_number] + log1p((u1 * g.rt_properties[0] - g.rt_properties[1]) * g.rt_properties[2]) * g.rt_properties[3];
        
        h_upper = g.rt_properties[4] * (sample - g.x[g.steps_number]) + g.rt_properties[5];
        
        u = unif_rand();
        
        if (u < f(sample, Rf, Renv) / exp(h_upper))
        {
          results[i] = sample;
          i++;
        }
        
        u1 = unif_rand();
        u1 = csl + u1 * (csr-csl);
        

      }else{
  
  u1 = (u1 - g.sampling_probabilities[0]) * g.unif_scaler;
  
  u1 *= g.steps_number;
  
  j = (int)u1;
  
  u1 -= j;
  
  if (u1 < g.p_a[j])
  { 
    u1 = u1 * g.s_upper_lower[j];
    
    sample = g.x[j] + u1 * (g.x[j + 1] - g.x[j]);
    
    results[i] = sample;
    
    i++;
    
    if (i < sample_size)
    {
      u1 = unif_rand();
      
      u1 = csl + u1 * (csr-csl);
    }
  }
  else
  {
    
    double u0 = unif_rand();
    
    sample = g.x[j] + u0 * (g.x[j + 1] - g.x[j]);
    
    f_sample = f(sample, Rf, Renv);
    
    double uf = f_sample /g.s_upper[j];
    
    if (u1 < uf)
    {
      
      results[i] = sample;
      i++;
    }
    
    u1 = unif_rand();
    u1 = csl + u1 * (csr-csl);
  }
  
}
      
}
  
  if(il != -1){
    g.x[il] = tmpxl;
  }

  
  if(ir != -1){
    g.x[ir] = tempxr;
  }

  
  PutRNGstate();
  
  UNPROTECT(1);
  
  return (Rresults);
  
}

