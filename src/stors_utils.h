#include "R_stors.h"

// return pdf value for the target dist
double f(double x, SEXP Rf, SEXP Renv)
{
  SEXP fcall, result, foo, val;
  double res;
  
  val = PROTECT(allocVector(REALSXP, 1));
  
  REAL(val)
    [0] = x;
  
  PROTECT(fcall = Rf_lang2(Rf, val));
  PROTECT(result = eval(fcall, Renv));
  PROTECT(foo = coerceVector(result, REALSXP));
  res = REAL(foo)[0];
  UNPROTECT(4);
  return (res);
}








