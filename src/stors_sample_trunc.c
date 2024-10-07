
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#if defined(CNUM) && defined(NAME)

#define IN_LEFT_TAIL -1
#define IN_RIGHT_TAIL -1


SEXP DEN_TRUNC(NAME)(SEXP s_size, SEXP Rxl, SEXP Rxr , SEXP Rcsl, SEXP Rcsr,  SEXP Ril, SEXP Rir){

    
    int j, sample_size = asInteger(s_size), il = asInteger(Ril), ir = asInteger(Rir);
    double csl = asReal(Rcsl), csr = asReal(Rcsr), tempxl, tempxr;
    double xl = asReal(Rxl), xr = asReal(Rxr);
    
    if(grids.grid[CNUM].x == NULL){
      REprintf("you need to optimize your destribution grid first");
      R_RETURN_NULL
    }
      
    struct grid g = grids.grid[CNUM];
      
#ifdef L_TAIL
    int check_lower = 0;

    if(il != IN_LEFT_TAIL){
      tempxl = g.x[il];
        g.x[il] = xl;
    }else{
      if(csl == 0) check_lower = 1;
    }
    
#else
    tempxl = g.x[il];
    g.x[il] = xl;

#endif

#ifdef R_TAIL
    int check_upper = 0;

      if(ir != IN_RIGHT_TAIL){
        tempxr = g.x[ir];
          g.x[ir] = xr;
          }else{
            if(csr == 1) check_upper = 1;
          }
#else
          tempxr = g.x[ir];
          g.x[ir] = xr;
#endif

      
#if L_TAIL == ARS || R_TAIL == ARS 
      
      double h_upper, u;
      
#endif
      
      double  u1, sample, f_sample;
      
      SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));
      
      double *results = REAL(Rresults);
      
      GetRNGstate();
      
      u1 = unif_rand();
      u1 = csl + u1 * (csr-csl);

      for (int i = 0; i < sample_size;)
      {
        
        
#ifdef L_TAIL    
        
        if (u1 < g.sampling_probabilities[0])
        {
          
#if L_TAIL == IT
          
          results[i] = L_ITF(u1);
          
          if(check_lower){
            if(results[i] >=  xl) i++;
          }else{
            i++;
          }
          
          u1 = unif_rand();
          u1 = csl + u1 * (csr-csl);

#elif L_TAIL == ARS
          
          sample = g.x[0] + (log( g.lt_properties[0] + u1 * g.lt_properties[1]) - g.lt_properties[2]) * g.lt_properties[3];
          h_upper = g.lt_properties[4] * (sample - g.x[0]) + g.lt_properties[2];
          u = unif_rand();
          
          if (u < F(sample) / exp(h_upper))
          {
            results[i] = sample;
            
            if(check_lower){
              if(results[i] >=  xl) i++;
            }else{
              i++;
            }
          }
          
          u1 = unif_rand();
          
          u1 = csl + u1 * (csr-csl);

#endif
          

          
        }else
          
#endif
          
          
#ifdef R_TAIL    
          
          if(u1 > g.sampling_probabilities[1]){
            
#if R_TAIL == IT
            
            results[i] = R_ITF(u1);
            
            if(check_upper){
              if(results[i] <=  xr) i++;
            }else{
              i++;
            }
            
            u1 = unif_rand();
            u1 = csl + u1 * (csr-csl);
            


#elif R_TAIL == ARS
            
            sample = g.x[g.steps_number] + log1p((u1 * g.rt_properties[0] - g.rt_properties[1]) * g.rt_properties[2]) * g.rt_properties[3];
            
            h_upper = g.rt_properties[4] * (sample - g.x[g.steps_number]) + g.rt_properties[5];
            
            u = unif_rand();
            
            if (u < F(sample) / exp(h_upper))
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
    
    f_sample = F(sample);
    
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
      
      
      if(il != IN_LEFT_TAIL){
        g.x[il] = tempxl;
      }

      if(ir != IN_RIGHT_TAIL){
        g.x[ir] = tempxr;
      }

      
      PutRNGstate();
      
      UNPROTECT(1);
      
      return (Rresults);
      
  }
  
  
#endif
