
#include "stors.h"

#include "cache.h"

#include "stors_utils.h"




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
  int check_lower = 0, check_upper =0;
  int Cnum = asInteger(Rcnum);
  
  
  struct grid g = grids.grid[Cnum];
  
  if(g.exist == 0){
    REprintf("you need to optimize your destribution grid first");
    R_RETURN_NULL
  }
  
  
  if(il != -1){
    tmpxl = g.x[il];
    g.x[il] = xl;
    if(csl == 0) check_lower = 1;
  }
  
  if(ir != -1){
    tempxr =g.x[ir];
    g.x[ir] = xr;
    if(csr == 1) check_upper = 1;
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
      
        if( check_lower ) {
          if(results[i] >=  xl) i++;
        }else{
          i++;
        }
        
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
        
        if(check_upper){
          if(results[i] <=  xr) i++;
        }else{
          i++;
        }
        
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

