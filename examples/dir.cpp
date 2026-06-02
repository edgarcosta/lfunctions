/*
Make up a degree 2 L-function by multiplying one of the non-real primitive
characters mod 5 and the quadratic character mod 7 together.
*/

#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file (as ./build/examples/dir.exe)
 * should generate something comparable to:
Command Line:- ./build/examples/dir.exe
Order of vanishing = 0
Sign = (0.85065080835203993218 + 0.52573111211913360603j)  +/-  (2.99e-58, 5.67e-58j)
First non-zero Taylor coeff = 0.91035185922722505047 +/- 4.7951e-57
L(1) = (1.0268799643452569392 + 0.24241347631804096081j)  +/-  (1.17e-38, 1.17e-38j)
First 10 zeros
Zero 0 = 4.4757382837286831320 +/- 2.0881e-53
Zero 1 = 6.1835781954508539144 +/- 6.6819e-52
Zero 2 = 6.8454917124913772678 +/- 1.3364e-51
Zero 3 = 8.4572291744232307216 +/- 2.6728e-51
Zero 4 = 11.160184543119529655 +/- 1.7106e-49
Zero 5 = 12.489603343033134238 +/- 5.4738e-48
Zero 6 = 12.674946417011355780 +/- 1.0948e-47
Zero 7 = 14.825025570328428251 +/- 3.5032e-46
Zero 8 = 15.112882258743768660 +/- 7.0065e-46
Zero 9 = 16.802876475728853982 +/- 2.8026e-45
Z-plot in [0, 10]:
0.00	0.91	                              |-----o
0.50	1.05	                              |------o
1.00	1.38	                              |---------o
1.50	1.87	                              |-------------o
2.00	2.40	                              |-----------------o
2.50	2.79	                              |-------------------o
3.00	2.78	                              |-------------------o
3.50	2.21	                              |---------------o
4.00	1.14	                              |-------o
4.48	zero	                              Z
4.50	-0.05	                              o
5.00	-0.83	                        o-----|
5.50	-0.82	                        o-----|
6.00	-0.22	                             o|
6.18	zero	                              Z
6.50	0.20	                              |o
6.85	zero	                              Z
7.00	-0.24	                             o|
7.50	-1.31	                     o--------|
8.00	-1.64	                  o-----------|
8.46	zero	                              Z
8.50	0.27	                              |o
9.00	4.33	                              |------------------------------
9.50	8.22	                              |------------------------------
10.00	8.97	                              |------------------------------
 */

#define __STDC_FORMAT_MACROS
#include <chrono>
#include <cstdint>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <flint/fmpz.h>
#include <flint/acb_poly.h>
#include <cassert>
#include "glfunc.h"


using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;

// compute the Euler poly for p
// with L the product of non-principal characters mod 5 and 7
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  acb_poly_t p5;
  acb_poly_init(p5);
  acb_poly_one(p5);

  /*
  if((p%5==1)||(p%5==4))
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5==2)||(p%5==3))
    acb_poly_set_coeff_si(p5,1,1);
  */

  // the Euler polynomials are (1-chi(p)p^-s)^-1
  acb_t i;
  acb_init(i);
  arb_set_ui(acb_imagref(i),1); // i
  if((p%5)==1)
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5)==2)
    acb_poly_set_coeff_acb(p5,1,i);
  if((p%5)==3)
    {
      acb_neg(i,i);
      acb_poly_set_coeff_acb(p5,1,i);
      acb_neg(i,i);
    }
  if((p%5)==4)
    acb_poly_set_coeff_si(p5,1,1);
  acb_clear(i);

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
  double mus[]={1,1}; // both dirichlet characters are odd
  Lerror_t ecode;

  // we have a degree 2 L-function with cond=5*7, alg=anal so
  // normalisation = 0.0
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
  printf("Order of vanishing = %" PRIu64 "\n",Lfunc_rank(L));
  printf("Sign = ");
  acb_printd(Lfunc_sign(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_sign(L) << endl;
  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_Taylor(L) << endl;


  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 1.0, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  { // regression assert: L(1) (complex)
    acb_t ref; acb_init(ref);
    arb_set_str(acb_realref(ref), "1.0268799643452569392", 300);
    arb_set_str(acb_imagref(ref), "0.24241347631804096081", 300);
    arb_add_error_2exp_si(acb_realref(ref), -50);
    arb_add_error_2exp_si(acb_imagref(ref), -50);
    assert(acb_overlaps(ctmp, ref));
    acb_clear(ref);
  }
  ecode|=Lfunc_special_value(ctmp, L, 2,0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  acb_clear(ctmp);

  printf("First 10 zeros\n");
  // we could use Lfunc_zeros(L, 1) for the dual L-function
  arb_srcptr zeros=Lfunc_zeros(L, 0);
  for(int i  = 0; i < 10; ++i) {
    printf("Zero %d = ", i);
    arb_printd(zeros+i, DIGITS);
    printf("\n");
    if (RAW) cout<<"RAW: "<<zeros + i<< endl;
  }
  { // regression assert: first zero
    arb_t ref; arb_init(ref);
    arb_set_str(ref, "4.4757382837286831320", 300);
    arb_add_error_2exp_si(ref, -50);
    assert(arb_overlaps(zeros + 0, ref));
    arb_clear(ref);
  }

  printf("Z-plot in [0, 10]:\n");
  Lplot_t *Lpp=Lfunc_plot_data(L, 0, 10.0, 20);
  int z = 0;
  double zero_double = arf_get_d(arb_midref(zeros + z), ARF_RND_NEAR);
  for(size_t k=0; k < Lpp->n_points; ++k) {
    printf("%.2f\t%.2f\t", k*Lpp->spacing , Lpp->points[k]);
    int y = 30 + int(7.5*Lpp->points[k]);
    int zero = 30;
    // assuming 60 columns
    for(int i = 0; i < 61; ++i) {
      if(i == y) {
        printf("o");
      } else if (i == zero) {
        printf("|");
      } else if ( (i > zero and i < y) or (i < zero and i > y) ) {
        printf("-");
      } else {
        printf(" ");
      }
    }
    printf("\n");
    if(k*Lpp->spacing < zero_double and (k+1)*Lpp->spacing >= zero_double){
      printf("%.2f\tzero\t", zero_double);
      for(int i = 0; i < 30; ++i)
        printf(" ");
      printf("Z\n");
      zero_double = arf_get_d(arb_midref(zeros + ++z), ARF_RND_NEAR);
    }
  }

  //free memory
  Lfunc_clear_plot(Lpp);
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}


