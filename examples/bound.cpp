#define __STDC_FORMAT_MACROS
#include "assert.h"
#include "acb_poly.h"
#include "glfunc.h"
#include "glfunc_internals.h"



int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i = 1; i < argc; i++)
    printf(" %s",argv[i]);
  printf("\n");

  uint64_t d = atoi(argv[1]);
  assert( d > 0 );
  assert((uint64_t)argc == d + 3);
  double normalisation = atof(argv[2]);
  printf("%" PRIu64 ":%.2f:[", d, normalisation);
  double * mus = new double[d];
  for(size_t i = 0; i < d; ++i) {
    mus[i] = atof(argv[3 + i]);
    printf("%f", mus[i]);
    if(i < d - 1) {
      printf(",");
    } else {
      printf("]");
    }
  }


  Lfunc_t Lf;
  Lerror_t ecode;

  //normalisation doesn't play a role
  Lf = Lfunc_init(d, 1, normalisation, mus, &ecode);
  if( fatal_error(ecode) ) {
    fprint_errors(stderr,ecode);
    return 1;
  }
  __attribute__((unused)) uint64_t target_M = Lfunc_nmax(Lf) ; // how many Euler polys does program want

  Lfunc *L;
  L=(Lfunc *)Lf;

  printf(":%.2f\n", ceil(100*exp(2*M_PI*(L->hi_i+0.5)*L->one_over_B))/100);

  Lfunc_clear(Lf);

  // print any warnings collected along the way
  fprint_errors(stderr, ecode);

  delete[] mus;
  return 0;
}


