
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"


#if defined(CNUM) && defined(NAME)

SEXP DEN_SAMPLE(NAME)(SEXP s_size, SEXP Rpassed_params){

  int j, sample_size = asInteger(s_size),
    match = TRUE, scale_at_end = FALSE;

  if(grids.grid[CNUM].x == NULL){
    REprintf("you need to optimize your destribution's grid first");
    R_RETURN_NULL
  }
  
  int n_param = grids.grid[CNUM].n_params;
  
  struct grid g = grids.grid[CNUM];

#ifdef SCALABLE

  double *pp = REAL(Rpassed_params);
  
  for(int i = 0; i < n_param; i++){
    if(g.params[i] != pp[i]) {match = FALSE; break;}
  }
  
  if(match == FALSE){
    
    match = TRUE;
    
    if(grids.grid[CNUM_SCALED].x != NULL){
      
      g = grids.grid[CNUM_SCALED];
      
      for(int i = 0; i < n_param; i++){
        if(g.params[i] != pp[i]) {match = FALSE; break;}
      }
      
      if(match == FALSE){
        scale_at_end = TRUE;
        g = grids.grid[CNUM];
      }
      
    }else{
      
      scale_at_end = TRUE;
      g = grids.grid[CNUM];
    }

  }

#endif
  

  
  
#if L_TAIL == ARS || R_TAIL == ARS 
  
  double h_upper, u;
 
#endif

  int flip = 9999;
  
  double  u1, sample, f_sample;
  
  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));
  
  double *results = REAL(Rresults);
  
  
#ifdef SPECIAL_FUNCTION
  SPECIAL_FUNCTION(sample_size)
#endif
  
  GetRNGstate();
  
  u1 = unif_rand();
  
  for (int i = 0; i < sample_size;)
  {
    

#ifdef FLIP_SAMPLE
    if(g.is_symmetric == 1){
      if(u1 > 0.5){
         u1 = 1-u1;
        flip = TRUE; 
      }else{
        flip = FALSE;
      }
    }
#endif
    

#ifdef L_TAIL    
    
    if (u1 < g.sampling_probabilities[0])
    {
      
#if L_TAIL == IT
      
#ifdef FLIP_SAMPLE
      if(g.is_symmetric == 1){
        results[i] = FLIP_SAMPLE(L_ITF(u1) ,flip);
        
      }else{
        results[i] = L_ITF(u1);    
      }
      
#else
      results[i] = L_ITF(u1);
      
#endif

      i++;
      u1 = unif_rand();
      
      
#elif L_TAIL == ARS
      
      sample = g.x[0] + (log( g.lt_properties[0] + u1 * g.lt_properties[1]) - g.lt_properties[2]) * g.lt_properties[3];
      h_upper = g.lt_properties[4] * (sample - g.x[0]) + g.lt_properties[2];
      u = unif_rand();
      if (u < F(sample) / exp(h_upper))
      {

#ifdef FLIP_SAMPLE
        if(g.is_symmetric == 1){
          results[i] = FLIP_SAMPLE(sample, flip);
        }else{
          results[i] = sample;    
        }
#else
        results[i] = sample;
#endif
        
        i++;
      }
      
      u1 = unif_rand();
      
      
#endif
      
    }else
      
#endif
      

#ifdef R_TAIL    

if(u1 > g.sampling_probabilities[1]){
      
#if R_TAIL == IT
      
      results[i] = R_ITF(u1);
      i++;
      u1 = unif_rand();
      

#elif R_TAIL == ARS
      
      sample = g.x[g.steps_number] + log1p((u1 * g.rt_properties[0] - g.rt_properties[1]) * g.rt_properties[2]) * g.rt_properties[3];
      
      h_upper = g.rt_properties[4] * (sample - g.x[g.steps_number]) + g.rt_properties[5];
      
      u = unif_rand();
      
      if (u < F(sample) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }
      
      u1 = unif_rand();
      
      
#endif
      
    }else
      
#endif

      {
      
      u1 = (u1 - g.sampling_probabilities[0]) * g.unif_scaler;
      
      u1 *= g.steps_number;
      
      j = (int)u1;
      
      u1 -= j;
      
      if (u1 < g.p_a[j])
      { 
        u1 = u1 * g.s_upper_lower[j];
        
        sample = g.x[j] + u1 * (g.x[j + 1] - g.x[j]);
        
#ifdef FLIP_SAMPLE

        if(g.is_symmetric == 1){
          results[i] = FLIP_SAMPLE(sample,flip);
        }else{
          results[i] = sample;    
        }
#else
        results[i] = sample;    
#endif
        
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
        
        f_sample = F(sample);
        
        double uf = f_sample /g.s_upper[j];
        
        if (u1 < uf)
        {
          
#ifdef FLIP_SAMPLE

          if(g.is_symmetric == 1){
            results[i] = FLIP_SAMPLE(sample,flip);
          }else{
            results[i] = sample;    
          }
#else
          results[i] = sample;    
#endif
          
          i++;
        }
        
        u1 = unif_rand();
        
      }
      
    }
    
  }
  
  
#ifdef SCALE
  if(scale_at_end){
    for(int i=0;i < sample_size; i++ ){
      results[i] = SCALE(results[i]);
    }
  }
#endif
  
  PutRNGstate();
  
  UNPROTECT(1);
  
  return (Rresults);
  
}


#endif
