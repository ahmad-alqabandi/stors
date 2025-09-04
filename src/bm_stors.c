#include "bm_stors.h"

#include "cache.h"

#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "uniform_rngs.h"

SEXP stors_srnorm_custom_inplace(SEXP Rresults){

  struct grid *g = grids.grid + CNUM + 1;
  double * restrict p_a = g->p_a;
  double * restrict x = g->x;

  int sample_size = LENGTH(Rresults);

  double *results = REAL(Rresults);

  int j;

  double h_upper, u;


  double  u1, sample, f_sample;

  GetRNGstate();

  u1 = u_rng();

  for (int i = 0; i < sample_size;)
  {



    if (u1 < g->sampling_probabilities[0])
    {

      sample = x[0] + (log( g->lt_properties[0] + u1 * g->lt_properties[1]) - g->lt_properties[2]) * g->lt_properties[3];
      h_upper = g->lt_properties[4] * (sample - x[0]) + g->lt_properties[2];
      u = u_rng();

      if (u < F(sample) / exp(h_upper))
      {
        results[i] = sample;
        i++;
      }

      u1 = u_rng();



    }else if(u1 > g->sampling_probabilities[1]){



        sample = x[g->steps_number] + log1p((u1 * g->rt_properties[0] - g->rt_properties[1]) * g->rt_properties[2]) * g->rt_properties[3];

        h_upper = g->rt_properties[4] * (sample - x[g->steps_number]) + g->rt_properties[5];

        u = u_rng();

        if (u < F(sample) / exp(h_upper))
        {
          results[i] = sample;
          i++;
        }

        u1 = u_rng();


      }else{

        u1 = (u1 - g->sampling_probabilities[0]) * g->unif_scaler;

        u1 *= g->steps_number;

        j = (int) u1;

        u1 -= j;

        if (u1 < p_a[j])
        {
          u1 = u1 * g->s_upper_lower[j];

          sample = x[j] + u1 * (x[j + 1] - x[j]);

          results[i] = sample;
          i++;

          if (i < sample_size)
          {
            u1 = u_rng();
          }

        }else{

          double u0 = u_rng();

          sample = x[j] + u0 * (x[j + 1] - x[j]);

          f_sample = F(sample);

          double uf = f_sample /g->s_upper[j];

          if (u1 < uf)
          {

            results[i] = sample;
            i++;
          }

          u1 = u_rng();

        }

      }

  }


  PutRNGstate();


  return (Rresults);


}