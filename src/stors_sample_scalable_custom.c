
#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#if defined(CUSTOM) || defined(SCALABLE) && defined(NAME)

// FOR NON-SYMMETRIC DIST

#ifndef FLIP_SAMPLE

#ifdef SCALABLE
#ifndef INPLACE
SEXP DEN_SAMPLE_SCALED(NAME)(SEXP s_size, SEXP Rpassed_params){

  struct grid *restrict g = grids.grid + CNUM ;
  double * restrict p_a = g->p_a;
  double * restrict x = g->x;
  int sample_size = asInteger(s_size);

  int match = TRUE;
  double *pp = REAL(Rpassed_params);
  int n_params = g->n_params;

  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));

  double *results = REAL(Rresults);
#endif

        #ifdef INPLACE
        SEXP DEN_SAMPLE_SCALED_INPLACE(NAME)( SEXP Rpassed_params, SEXP Rresults){

          struct grid *restrict g = grids.grid + CNUM ;
          double * restrict p_a = g->p_a;
          double * restrict x = g->x;

          int match = TRUE;
          double *pp = REAL(Rpassed_params);
          int n_params = g->n_params;

          int sample_size = LENGTH(Rresults);

          double *results = REAL(Rresults);

  #endif
#endif



#ifdef CUSTOM
#ifndef INPLACE

  SEXP DEN_SAMPLE_CUSTOM(NAME)(SEXP s_size){
    struct grid *g = grids.grid + CNUM + 1;
    double * restrict p_a = g->p_a;
    double * restrict x = g->x;
    int sample_size = asInteger(s_size);

    SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));

    double *results = REAL(Rresults);

#endif

    #ifdef INPLACE
    SEXP DEN_SAMPLE_CUSTOM_INPLACE(NAME)(SEXP Rresults){
      struct grid *g = grids.grid + CNUM + 1;
      double * restrict p_a = g->p_a;
      double * restrict x = g->x;

      int sample_size = LENGTH(Rresults);

      double *results = REAL(Rresults);

  #endif
#endif
#endif


// FOR SYMMETRIC DIST
#ifdef FLIP_SAMPLE

#ifdef SCALABLE
#ifndef INPLACE
SEXP DEN_SAMPLE_SYM_SCALED(NAME)(SEXP s_size, SEXP Rpassed_params){

  struct grid *restrict g = grids.grid + CNUM ;
  double * restrict p_a = g->p_a;
  double * restrict x = g->x;
  int sample_size = asInteger(s_size);

  int match = TRUE;
  double *pp = REAL(Rpassed_params);
  int n_params = g->n_params;

  SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));

  double *results = REAL(Rresults);
#endif

  #ifdef INPLACE
        SEXP DEN_SAMPLE_SYM_SCALED_INPLACE(NAME)( SEXP Rpassed_params, SEXP Rresults){

          struct grid *restrict g = grids.grid + CNUM ;
          double * restrict p_a = g->p_a;
          double * restrict x = g->x;

          int match = TRUE;
          double *pp = REAL(Rpassed_params);
          int n_params = g->n_params;

          int sample_size = LENGTH(Rresults);

          double *results = REAL(Rresults);

  #endif

#endif

#ifdef CUSTOM
#ifndef INPLACE

  SEXP DEN_SAMPLE_SYM_CUSTOM(NAME)(SEXP s_size){
    struct grid *g = grids.grid + CNUM + 1;
    double * restrict p_a = g->p_a;
    double * restrict x = g->x;
    int sample_size = asInteger(s_size);

    SEXP Rresults = PROTECT((allocVector(REALSXP, sample_size)));

    double *results = REAL(Rresults);
#endif
  #ifdef INPLACE
    SEXP DEN_SAMPLE_SYM_CUSTOM_INPLACE(NAME)(SEXP Rresults){
      struct grid *g = grids.grid + CNUM + 1;
      double * restrict p_a = g->p_a;
      double * restrict x = g->x;

      int sample_size = LENGTH(Rresults);

      double *results = REAL(Rresults);

    #endif
  #endif
#endif

  int j;

#if L_TAIL == ARS || R_TAIL == ARS

  double h_upper, u;

#endif

#ifdef FLIP_SAMPLE

  int flip;

#endif


  double  u1, sample, f_sample;


#ifdef SPECIAL_FUNCTION
  SPECIAL_FUNCTION(sample_size)
#endif

  GetRNGstate();

  u1 = unif_rand();

  for (int i = 0; i < sample_size;)
  {


#ifdef FLIP_SAMPLE
      if(u1 > 0.5){
         u1 = 1-u1;
        flip = TRUE;
      }else{
        flip = FALSE;
      }
#endif


#ifdef L_TAIL

    if (u1 < g->sampling_probabilities[0])
    {

#if L_TAIL == IT

#ifdef FLIP_SAMPLE
        results[i] = FLIP_SAMPLE(L_ITF(u1 * g->proposal_area + CDF(g->lower)) ,flip);
#else
      double U = u1 * g->proposal_area + CDF(g->lower);
      results[i] = L_ITF(U);

#endif

      i++;
      u1 = unif_rand();


#elif L_TAIL == ARS

      sample = x[0] + (log( g->lt_properties[0] + u1 * g->lt_properties[1]) - g->lt_properties[2]) * g->lt_properties[3];
      h_upper = g->lt_properties[4] * (sample - x[0]) + g->lt_properties[2];
      u = unif_rand();

      if (u < F(sample) / exp(h_upper))
      {

#ifdef FLIP_SAMPLE
          results[i] = FLIP_SAMPLE(sample, flip);
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

if(u1 > g->sampling_probabilities[1]){

#if R_TAIL == IT

      double U = CDF(g->upper) - g->proposal_area + u1 * g->proposal_area;
      results[i] = R_ITF(U);

      i++;
      u1 = unif_rand();

#elif R_TAIL == ARS

      sample = x[g->steps_number] + log1p((u1 * g->rt_properties[0] - g->rt_properties[1]) * g->rt_properties[2]) * g->rt_properties[3];

      h_upper = g->rt_properties[4] * (sample - x[g->steps_number]) + g->rt_properties[5];

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

      u1 = (u1 - g->sampling_probabilities[0]) * g->unif_scaler;

      u1 *= g->steps_number;

      j = (int) u1;

      u1 -= j;

      if (u1 < p_a[j])
      {
        u1 = u1 * g->s_upper_lower[j];

        sample = x[j] + u1 * (x[j + 1] - x[j]);

#ifdef FLIP_SAMPLE
        results[i] = FLIP_SAMPLE(sample,flip);
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

        sample = x[j] + u0 * (x[j + 1] - x[j]);

        f_sample = F(sample);

        double uf = f_sample /g->s_upper[j];

        if (u1 < uf)
        {

#ifdef FLIP_SAMPLE

          results[i] = FLIP_SAMPLE(sample,flip);
#else
          results[i] = sample;
#endif

          i++;
        }

        u1 = unif_rand();

      }

    }

  }


#ifdef SCALABLE

  for(int i=0; i<n_params; i++){
    if(g->params[i] != pp[i]){
      match = FALSE;
      break;
    }
  }

  if(!match){
    for(int i=0;i < sample_size; i++ ){
      results[i] = SCALE(results[i]);
    }
  }

#endif

  PutRNGstate();

#ifndef INPLACE
  UNPROTECT(1);
#endif

  return (Rresults);

}


#endif
