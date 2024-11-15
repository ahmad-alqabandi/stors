
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#if defined(CNUM) && defined(NAME)

SEXP DEN_TRUNC_NAV(NAME)(SEXP Rlx, SEXP Rrx){
  
if(grids.grid[CNUM].x == NULL){
  REprintf("you need to optimize your destribution grid first");
  R_RETURN_NULL
}

struct grid g = grids.grid[CNUM];

int i=0;

double xlr[2], total_area = g.areas[0] + g.areas[1]+ g.areas[2], ixlr[2];


#if R_TAIL == ARS

double hu_x;

#elif R_TAIL == IT


#endif


xlr[0] = asReal(Rlx);
xlr[1] = asReal(Rrx);

SEXP Rresults = PROTECT((allocVector(REALSXP, 4)));

double *results = REAL(Rresults);

for( int j = 0; j < 2; j++){ 
  
  
  
#ifdef R_TAIL 
  
  if(xlr[j] > g.x[g.steps_number]){
    
#if R_TAIL == ARS
    
    hu_x =  g.rt_properties[4] * (xlr[j] - g.x[g.steps_number]) + g.rt_properties[5];
    
    results[j] =( g.areas[0] + g.areas[1] + (g.rt_properties[3] * (exp(hu_x) - exp(g.rt_properties[5]))) )/ total_area;
    
    ixlr[j] = -1;
    
#elif R_TAIL == IT

    
    //cdf = CDF(xlr[j]) - CDF(g.x[g.steps_number]);

    //results[j] =  ( g.areas[0] + g.areas[1] + cdf)/ total_area;
    
    // added line _ this is because in tails we use IT ( the target IT function)
    results[j]  = CDF(xlr[j]); 
    
    ixlr[j] = -1;
    
#endif
    
  }else
    
#endif

    {
    while( i < g.steps_number + 1){
      
      if(xlr[j] < g.x[i]){
        
        
#ifdef L_TAIL
        
    if(i == 0){

#if L_TAIL == ARS

          hu_x =  g.lt_properties[4] * (xlr[j] - g.x[0]) + g.lt_properties[2];
          
          results[j] = (g.lt_properties[3] * (exp(hu_x) - g.lt_properties[0])) / total_area;
          ixlr[j] = -2;
#elif L_TAIL == IT
          
          results[j] = CDF(xlr[j]);
            
          ixlr[j] = -2;
#endif
          
        }else
          
#endif
    
          {
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
#endif
