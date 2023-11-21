#include "R_stors.h"
#include "R_cache.h"


struct grids grids = {{0}};



SEXP cache_grid(SEXP R_Cnum, SEXP R_x, SEXP R_s_upper, SEXP R_p_a, SEXP R_s_upper_lower, SEXP R_areas, SEXP R_steps_number, SEXP R_sampling_probabilities, SEXP R_unif_scaler, SEXP R_lt_properties, SEXP R_rt_properties){
  
  int j = asInteger(R_Cnum), m = asInteger(R_steps_number);
  double unif_scaler = asReal(R_unif_scaler);
  double *x = REAL(R_x), *S_upper = REAL(R_s_upper), *p_a = REAL(R_p_a), *s_upper_lower = REAL(R_s_upper_lower),
    *areas = REAL(R_areas), *sampling_probabilities = REAL(R_sampling_probabilities), *lt_properties = REAL(R_lt_properties), *rt_properties = REAL(R_rt_properties);
  
  grids.grid[j].x =  calloc( (m + 1) , sizeof(double) );
  grids.grid[j].p_a = calloc( (m + 1) , sizeof(double) );
  grids.grid[j].s_upper = calloc( (m + 1) , sizeof(double) );
  grids.grid[j].s_upper_lower = calloc( (m + 1) ,  sizeof(double) );
  
  grids.grid[j].steps_number = m;
  grids.grid[j].unif_scaler = unif_scaler;

  for( size_t i = 0; i < (m + 1); i++){
    grids.grid[j].x[i] = x[i];
    grids.grid[j].p_a[i] = p_a[i];
    grids.grid[j].s_upper[i] = S_upper[i];
    grids.grid[j].s_upper_lower[i] = s_upper_lower[i];
  }
  
  for( size_t i = 0; i < 3; i++){
    grids.grid[j].areas[i] = areas[i];
  }
  
  for( size_t i = 0; i < 2; i++){
    grids.grid[j].sampling_probabilities[i] = sampling_probabilities[i];
  }
  
  for( size_t i = 0; i < 2; i++){
    grids.grid[j].sampling_probabilities[i] = sampling_probabilities[i];
  }
  
  for( size_t i = 0; i < 5; i++){
    grids.grid[j].lt_properties[i] = lt_properties[i];
  }
  
  for( size_t i = 0; i < 6; i++){
    grids.grid[j].rt_properties[i] = rt_properties[i];
  }
  
  grids.incache +=1;
  

  R_RETURN_NULL
}



SEXP free_cache(void){
  
  for( size_t i = 0; i < grids.incache; i++){
    free(grids.grid[i].x);
    free(grids.grid[i].p_a);
    free(grids.grid[i].s_upper);
    free(grids.grid[i].s_upper_lower);
  }
  
  grids.incache=0;
  
  Rprintf("\n === C Cache freed successfully === \n");
  
  R_RETURN_NULL
}

