#define __STDC_FORMAT_MACROS
#include "assert.h"
#include "acb_poly.h"
#include "glfunc.h"
#include "glfunc_internals.h"

int message(int argc, char **argv) {
  fprintf(stderr, "Usage:\n %s <degree> <normalisation> <mu_0> <mu_1> ... <mu_degree>\n with <degree> between 1 and 10\n", argv[0]);
  fprintf(stderr, "For example:\n");
  fprintf(stderr, "\t- L-function of an Elliptic curve over Q:\n");
  fprintf(stderr, "\t\t# %s 2 0.5 0 1\n", argv[0]);
  fprintf(stderr, "\t or\t# %s 2 0 0.5 1.5\n", argv[0]);
  fprintf(stderr, "\t- L-function of a classical modular form of weight w:\n");
  fprintf(stderr, "\t\t# %s 2 (w-1)*0.5 0 1\n", argv[0]);
  fprintf(stderr, "\t or\t# %s 2 0 (w-1)*0.5 (w-1)*0.5 + 1\n", argv[0]);
  if(argc >= 2) {
    fprintf(stderr, "Arguments given:\n");
    fprintf(stderr, "\t<degree> = %s\n", argv[1]);
  }
  if(argc >= 3)
    fprintf(stderr, "\t<normalisation> = %s\n", argv[2]);

  for(int i = 0; i < argc - 3; ++i)
    fprintf(stderr, "\t<mu_%d> = %s\n", i, argv[2]);
  return false;
}

int parse_arg(uint64_t &d, double &normalisation, double *mus, int argc, char **argv) {
  if(argc < 5)
    return message(argc, argv);

  d = atoi(argv[1]);
  if( (uint64_t)argc != d + 3 )
    return message(argc, argv);

  normalisation = atof(argv[2]);
  printf("%" PRIu64 ":%.2f:[", d, normalisation);
  for(size_t i = 0; i < d; ++i) {
    mus[i] = atof(argv[3 + i]);
    printf("%f", mus[i]);
    if(i < d - 1) {
      printf(",");
    } else {
      printf("]");
    }
  }
  return true;
}

int main (int argc, char**argv)
{

  uint64_t d;
  double normalisation;
  double *mus = new double[MAX_DEGREE];
  if( not parse_arg(d, normalisation, mus, argc, argv) )
    return 1;

  Lfunc_t Lf;
  Lerror_t ecode;

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


