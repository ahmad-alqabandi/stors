#include "R_cache.h"
#include "cache.h"


struct grids grids = {{0}};


SEXP cache_grid(SEXP R_Cnum, SEXP R_x, SEXP R_s_upper,
                SEXP R_p_a, SEXP R_s_upper_lower,
                SEXP R_areas, SEXP R_steps_number,
                SEXP R_sampling_probabilities, SEXP R_unif_scaler,
                SEXP R_lt_properties, SEXP R_rt_properties,
                SEXP Ralpha, SEXP Rsymmetric,
                SEXP Rparams, SEXP Rn_params, SEXP Rlb, SEXP Rrb){

  int j = asInteger(R_Cnum), m = asInteger(R_steps_number), n_params = asInteger(Rn_params);

  double unif_scaler = asReal(R_unif_scaler), alpha = asReal(Ralpha);

  double *x = REAL(R_x), *S_upper = REAL(R_s_upper), *p_a = REAL(R_p_a),
    *s_upper_lower = REAL(R_s_upper_lower), *areas = REAL(R_areas),
    *sampling_probabilities = REAL(R_sampling_probabilities),
    *lt_properties = REAL(R_lt_properties), *rt_properties = REAL(R_rt_properties),
    *params = REAL(Rparams), symmetric = 0;

    double lower = asReal(Rlb), upper = asReal(Rrb) ;

    if(!isNull(Rsymmetric)){
      symmetric = asReal(Rsymmetric);
      grids.grid[j].is_symmetric = 1;
    }else{
      grids.grid[j].is_symmetric = 0;
    }

  grids.grid[j].x =  R_Calloc( (m + 1) , double );
  grids.grid[j].p_a = R_Calloc( (m + 1) , double );
  grids.grid[j].s_upper = R_Calloc( (m + 1) , double );
  grids.grid[j].s_upper_lower = R_Calloc( (m + 1) ,  double );

  grids.grid[j].steps_number = m;
  grids.grid[j].unif_scaler = unif_scaler;
  grids.grid[j].alpha = alpha;
  grids.grid[j].symmetric = symmetric;
  grids.grid[j].n_params = n_params;
  grids.grid[j].lower = lower;
  grids.grid[j].upper = upper;
  grids.grid[j].proposal_area = 0;

  for( size_t i = 0; i < (m + 1); i++){
    grids.grid[j].x[i] = x[i];
    grids.grid[j].p_a[i] = p_a[i];
    grids.grid[j].s_upper[i] = S_upper[i];
    grids.grid[j].s_upper_lower[i] = s_upper_lower[i];
  }

  for( size_t i = 0; i < 3; i++){
    grids.grid[j].areas[i] = areas[i];
    grids.grid[j].proposal_area += areas[i];
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

  for( size_t i = 0; i < n_params; i++){
    grids.grid[j].params[i] = params[i];
  }

  grids.grid[j].exist = 1;

  grids.incache +=1;

  R_RETURN_NULL
}



SEXP free_cache(void){

  for( size_t i = 0; i < MAX_GRIDS_NUMBER; i++){

    if(grids.grid[i].x != NULL){

      R_Free(grids.grid[i].x);
      grids.grid[i].x = NULL;
      R_Free(grids.grid[i].p_a);
      grids.grid[i].p_a = NULL;
      R_Free(grids.grid[i].s_upper);
      grids.grid[i].s_upper = NULL;
      R_Free(grids.grid[i].s_upper_lower);
      grids.grid[i].s_upper_lower = NULL;

    }

  }

  grids.incache=0;

  Rprintf("\n === C Cache freed successfully === \n");

  R_RETURN_NULL
}



SEXP free_cache_cnum( SEXP Rcnum){

  int cnum = asInteger(Rcnum);

  R_Free(grids.grid[cnum].x);
  grids.grid[cnum].x = NULL;
  R_Free(grids.grid[cnum].p_a);
  grids.grid[cnum].p_a = NULL;
  R_Free(grids.grid[cnum].s_upper);
  grids.grid[cnum].s_upper = NULL;
  R_Free(grids.grid[cnum].s_upper_lower);
  grids.grid[cnum].s_upper_lower = NULL;

  grids.incache-=1;

  R_RETURN_NULL
}
