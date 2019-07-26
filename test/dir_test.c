/*
   Same as dir_test1.c but use callback interface instead.
*/

#include <inttypes.h>
#include "acb_poly.h"
#include "glfunc.h"

// compute the Euler poly for p
// with L the product of non-principal characters mod 5 and 7
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  // pretend we run out of polynomials at p>100
  if(p>100) {
    acb_poly_zero(poly);
    return;
  }
  acb_poly_t p5;
  acb_poly_init(p5);
  acb_poly_one(p5);
  if((p%5==1)||(p%5==4))
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5==2)||(p%5==3))
    acb_poly_set_coeff_si(p5,1,1);
  acb_poly_t p7;
  acb_poly_init(p7);
  acb_poly_one(p7);
  if((p%7==1)||(p%7==2)||(p%7==4))
    acb_poly_set_coeff_si(p7,1,-1);
  if((p%7==3)||(p%7==5)||(p%7==6))
    acb_poly_set_coeff_si(p7,1,1);
  acb_poly_mul(poly,p5,p7,prec);
  acb_poly_clear(p5);
  acb_poly_clear(p7);
}


int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");

  Lfunc_t L;
  double mus[]={0,1};
  Lerror_t ecode;

  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  L=Lfunc_init(2,5*7,0.0,mus,&ecode);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr, ecode);
    return 0;
  }

  // do the computation
  ecode|=Lfunc_compute(L);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  // now extract some information
  printf("(Apparent) Rank = %" PRIu64 "\n",Lfunc_rank(L));
  printf("Epsilon = ");acb_printd(Lfunc_epsilon(L),20);printf("\n");
  printf("First non-zero Taylor coeff = ");arb_printd(Lfunc_Taylor(L),20);printf("\n");

  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp,L,1.0,0.0);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }
  printf("L(1) = ");acb_printd(ctmp,20);printf("\n");
  acb_clear(ctmp);

  arb_srcptr zeros=Lfunc_zeros(L,0);
  for(uint64_t z=0;!arb_is_zero(zeros+z);z++)
  {printf("Zero %" PRIu64 " = ",z);arb_printd(zeros+z,20);printf("\n");}

  zeros=Lfunc_zeros(L,1);
  for(uint64_t z=0;!arb_is_zero(zeros+z);z++)
  {printf("Zero %" PRIu64 " = ",z);arb_printd(zeros+z,20);printf("\n");}

  Lplot_t *Lpp=Lfunc_plot_data(L,0,10.0,500);
  for(uint64_t i=0;i<Lpp->n_points;i++)
    printf("%f %f\n",i*Lpp->spacing,Lpp->points[i]);

  //free memory
  Lfunc_clear_plot(Lpp);
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}

