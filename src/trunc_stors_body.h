if(grids.grid[CNUM].x == NULL){
  REprintf("you need to optimize your destribution grid first");
  R_RETURN_NULL
}

struct grid g = grids.grid[CNUM];

int i=0;

double xlr[2], total_area = g.areas[0]+ g.areas[1]+ g.areas[2];


#if R_TAIL == ARS

double hu_x;

#elif R_TAIL == IT

double cdf;
// SEXP fcall, result, foo, val;
// double res;

#endif


xlr[0] = asReal(Rlx);
xlr[1] = asReal(Rrx);

SEXP Rresults = PROTECT((allocVector(REALSXP, 2)));

double *results = REAL(Rresults);

for( int j = 0; j < 2; j++){
  
  if(xlr[j] > g.x[g.steps_number]){
    
#if R_TAIL == ARS
    
    hu_x =  g.rt_properties[4] * (xlr[j] - g.x[g.steps_number]) + g.rt_properties[5];
    
    results[j] =( g.areas[0] + g.areas[1] + (g.rt_properties[3] * (exp(hu_x) - exp(g.rt_properties[5]))) )/ total_area;
    
#elif R_TAIL == IT

    
    // val = PROTECT(allocVector(REALSXP, 2));
    // 
    // REAL(val)[0] = g.x[g.steps_number];
    // REAL(val)[1] = xlr[j];
    // 
    // PROTECT(fcall = Rf_lang2(Rf, val));
    // PROTECT(result = eval(fcall, Renv));
    // PROTECT(foo = coerceVector(result, REALSXP));
    // cdf[0] = REAL(foo)[0];
    // cdf[1] = REAL(foo)[1];
    // UNPROTECT(4);
    
    cdf = CDF(xlr[j]) - CDF(g.x[g.steps_number]);

results[j] =  ( g.areas[0] + g.areas[1] + cdf)/ total_area;
    
#endif
    
  }else{
    while( i < g.steps_number + 1){
      
      if(xlr[j] <= g.x[i]){
        
        if(i == 0){
          
#if R_TAIL == ARS

          hu_x =  g.lt_properties[4] * (xlr[j] - g.x[0]) + g.lt_properties[2];
          
          results[j] = (g.lt_properties[3] * (exp(hu_x) - g.lt_properties[0])) / total_area;
          
#elif R_TAIL == IT
          
          
          // val = PROTECT(allocVector(REALSXP, 1));
          // 
          // REAL(val)[0] =  xlr[j];
          // 
          // PROTECT(fcall = Rf_lang2(Rf, val));
          // PROTECT(result = eval(fcall, Renv));
          // PROTECT(foo = coerceVector(result, REALSXP));
          // cdf[0] = REAL(foo)[0];
          // UNPROTECT(4);
          // 
          cdf = CDF(xlr[j]);
          
          results[j] =  cdf / total_area;
          
#endif
          
        }else{
          results[j] =(g.areas[0] + g.alpha * (i - 1) + (xlr[j] - g.x[i-1]) * g.s_upper[i-1]) / total_area;
        }
        break;
      }
      
      i++;
      
    }
  }
  
  
}

UNPROTECT(1);

return (Rresults);
