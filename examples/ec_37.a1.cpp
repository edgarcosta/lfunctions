// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Make up a degree 2 L-function associated the Elliptic Curve with LMFDB label 37.a1
 * https://www.lmfdb.org/EllipticCurve/Q/37/a/1
 *
 * Change the following two lines, to modify the number of decimal digits printed, or to print raw format (not human friendly)
 */
#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file should generate something
 * comparable to:
Order of vanishing = 1
Epsilon = (-1 + 4.1376012004041556323e-62j)  +/-  (4.98e-58, 6.37e-56j)
First non-zero Taylor coeff = 0.30599977383405230182 +/- 1.0943e-52
L(1.5) = (0.18396547525832984973 - 5.9947846803720119316e-58j)  +/-  (6.14e-39, 6.14e-39j)
L(2.5) = (0.55179233807261090057 - 1.6284871398197777388e-47j)  +/-  (5.66e-41, 5.66e-41j)
First 10 zeros
Zero 0 = 5.0031700140066586953 +/- 2.0881e-53
Zero 1 = 6.8703912169544319485 +/- 4.2764e-50
Zero 2 = 8.0143308078728792234 +/- 1.3364e-51
Zero 3 = 9.9330983536053517143 +/- 8.5528e-50
Zero 4 = 10.775138162540800446 +/- 6.8423e-49
Zero 5 = 11.757324722849776347 +/- 1.3685e-48
Zero 6 = 12.958386413882845948 +/- 2.7369e-48
Zero 7 = 15.603857873204318865 +/- 2.8026e-45
Zero 8 = 16.192017416874481955 +/- 2.8026e-45
Zero 9 = 17.141693648014874803 +/- 5.6052e-45
Z-plot in [0, 10]:
0.00	0.00	                              o
0.50	0.18	                              |o
1.00	0.48	                              |--o
1.50	1.02	                              |------o
2.00	1.83	                              |------------o
2.50	2.82	                              |--------------------o
3.00	3.74	                              |---------------------------o
3.50	4.15	                              |------------------------------
4.00	3.65	                              |--------------------------o
4.50	2.12	                              |--------------o
5.00	0.01	                              o
5.00	zero	                              Z
5.50	-1.73	                  o-----------|
6.00	-2.17	              o---------------|
6.50	-1.16	                      o-------|
6.87	zero	                              Z
7.00	0.36	                              |-o
7.50	0.99	                              |------o
8.00	0.05	                              o
8.01	zero	                              Z
8.50	-1.61	                  o-----------|
9.00	-2.30	             o----------------|
9.50	-1.33	                     o--------|
9.93	zero	                              Z
10.00	0.16	                              |o
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
#include <flint/fmpzxx.h>
#include <acb_poly.h>
#include "glfunc.h"
#include "examples_tools.h"

using flint::fmpzxx;
using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;


// such dictionary can be obtained directly from the sidebar
// http://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/23/1/b/a/
map<int64_t, vector<int64_t>> euler_factors  = {
  {2, {1, 2, 2}},
  {3, {1, 3, 3}},
  {5, {1, 2, 5}},
  {7, {1, 1, 7}},
  {11, {1, 5, 11}},
  {13, {1, 2, 13}},
  {17, {1, 0, 17}},
  {19, {1, 0, 19}},
  {23, {1, -2, 23}},
  {29, {1, -6, 29}},
  {31, {1, 4, 31}},
  {37, {1, 1}},
  {41, {1, 9, 41}},
  {43, {1, -2, 43}},
  {47, {1, 9, 47}},
  {53, {1, -1, 53}},
  {59, {1, -8, 59}},
  {61, {1, 8, 61}},
  {67, {1, -8, 67}},
  {71, {1, -9, 71}},
  {73, {1, 1, 73}},
  {79, {1, -4, 79}},
  {83, {1, 15, 83}},
  {89, {1, -4, 89}},
  {97, {1, -4, 97}},
  {101, {1, -3, 101}},
  {103, {1, -18, 103}},
  {107, {1, 12, 107}},
  {109, {1, 16, 109}},
  {113, {1, 18, 113}},
  {127, {1, -1, 127}},
  {131, {1, 12, 131}},
  {137, {1, 6, 137}},
  {139, {1, -4, 139}},
};


// compute the Euler poly for p
// just uses the map above
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec __attribute__((unused)), void *param __attribute__((unused)))
{
  acb_poly_zero(poly);
  auto it = euler_factors.find(p);
  if( it != euler_factors.end() ) {
  for(size_t i = 0; i < it->second.size(); ++i)
    acb_poly_set_coeff_si(poly, i, it->second[i]);
  }
}




int main ()
{
  Lfunc_t L;
  double mus[2] = {0, 1};
  Lerror_t ecode;

  // we have a degree 2 L-function of motivic weight 1, so normalisation = 0.5
  L = Lfunc_init(2, 37, 0.5, mus, &ecode);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  // populate local factors
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
  printf("Epsilon = ");
  acb_printd(Lfunc_epsilon(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_epsilon(L) << endl;
  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_Taylor(L) << endl;
  arb_t bsd;
  char bsd_str[] = "0.305999773834052301820483683321676474452637774590771998534541832481016050469290169911495257337795897237898682879524967997997869651621709648704953228700";
  arb_init(bsd);
  arb_set_str(bsd, bsd_str, 400);
  assert(arb_overlaps(Lfunc_Taylor(L), bsd));
  arb_clear(bsd);


  //FIXME add asserts and these two lines to intro text
  // ~ 0.183965475258329849732118662920
  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 1.5, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  // ~ 0.551792338072610900567950084313
  ecode|=Lfunc_special_value(ctmp, L, 2.5,0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  printf("L(2.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
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
