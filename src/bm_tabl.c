#include "bm_tabl.h"

#include "cache.h"

#include "R_stors.h"

#include "macro_func.h"

#include "macro_var.h"

#include "cache.h"

#include "uniform_rngs.h"

SEXP tabl_srnorm_custom_inplace(SEXP Rresults){

    struct grid *g = grids.grid + CNUM + 1;
    double * restrict p_a = g->p_a;
    double * restrict x = g->x;
    double tabl_pa, q, h, d, tabl_steps = g->steps_number - 1;

    int sample_size = LENGTH(Rresults);
    int flip;

    double *results = REAL(Rresults);

    int j, accept;

    double h_upper, u;




  double  u1,u2, sample, f_sample;

    GetRNGstate();

  u1 = u_rng();


  for (int i = 0; i < sample_size;)
  {


    // if(u1 > 0.5){
    //   u1 = 1-u1;
    //   flip = -1;
    // }else{
    //   flip = 1;
    // }

    /* Use the top bit for the sign, rescale the remainder back to (0,1)   */
    if (u1 < 0.5) {                    /* left half → negative sign          */
    flip = -1;
      u1 = 2.0 * u1;                /* uniform on (0,1) again             */
    } else {                          /* right half → positive sign         */
    flip = +1;
      u1 = 2.0 * u1 - 1.0;          /*   (shifts interval [0.5,1) → (0,1) */
    }
//
//     if (u1 < g->sampling_probabilities[0])
//     {
//
//       sample = x[0] + (log( g->lt_properties[0] + u1 * g->lt_properties[1]) - g->lt_properties[2]) * g->lt_properties[3];
//       h_upper = g->lt_properties[4] * (sample - x[0]) + g->lt_properties[2];
//       u = u_rng();
//
//       if (u < F(sample) / exp(h_upper))
//       {
//         results[i] = sample;
//         i++;
//       }
//
//       u1 = u_rng();
//
//
//
//     }else
//
//


      if(u1 > g->sampling_probabilities[1]){



        sample = x[g->steps_number] + log1p((u1 * g->rt_properties[0] - g->rt_properties[1]) * g->rt_properties[2]) * g->rt_properties[3];

        h_upper = g->rt_properties[4] * (sample - x[g->steps_number]) + g->rt_properties[5];

        u = u_rng();

        if (u < F(sample) / exp(h_upper))
        {
          results[i] = sample * flip;
          i++;
        }

        u1 = u_rng();


      }else
        // ====================== TABL ALGORITHM START ======================
        {


  u1 = (u1 - g->sampling_probabilities[0]) * g->unif_scaler;

  u1 *= tabl_steps;

  j = (int) u1;

  u1 -= j;

  h = (x[j+2] - x[j+1]);
  d = (x[j+1] - x[j]);


  // tabl_pa = h / d;

  q = u1 * h;

  if (q <= d) {
    results[i++] = (x[j] + q) * flip;
    u1 = u_rng();
    continue;
  }

    double u0 = u_rng();

    sample = (x[j] + u0 * (x[j + 1] - x[j])) * flip;


    f_sample = F(sample);

    double uf = f_sample /g->s_upper[j];

    if (u1 < uf)
    {

      results[i] = sample;
      i++;
    }

    u1 = u_rng();


}
      // ====================== TABL ALGORITHM END ======================

  }


  PutRNGstate();


  return (Rresults);


}