#include "tstors.h"


#define add(x) (x+10)
#define sub(x) (x-10)

SEXP calc(SEXP Rnum){
  
  int num = asInteger(num);
  
  if( DEF == 1 ){
    Rprintf(" \n\n result = %d", add(num));
  }else{
    Rprintf(" \n\n result = %d", sub(num));
  }
  
  R_RETURN_NULL
}