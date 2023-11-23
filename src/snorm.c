#include "R_stors.h"
#include "snorm.h"

SEXP srnorm(SEXP s_size)
{
  int j, sample_size = asInteger(s_size);
  
   if(grids.grid[0].x == NULL){
     REprintf("you need to optimize your destribution grid first");
     R_RETURN_NULL
   }
  
  struct grid g = grids.grid[0];
  
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
      
      if (u < 0.3989423 * exp(-0.5 * sample * sample) / exp(h_upper))
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
      
      if (u < 0.3989423 * exp(-0.5 * sample * sample) / exp(h_upper))
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
        
        f_sample =  0.3989423 * exp(-0.5 * sample * sample);
        
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

