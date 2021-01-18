// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Make up a degree 2 L-function associated to the Ramanujan tau function.
 * See:
 * - https://en.wikipedia.org/wiki/Ramanujan_tau_function
 * - https://www.lmfdb.org/ModularForm/GL2/Q/holomorphic/1/12/a/a/
 *
 * Change the following two lines, to modify the number of decimal digits printed, or to print raw format (not human friendly)
 */
#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file should generate something
 * comparable to:
Order of vanishing = 0
Epsilon = (1 + 0j)  +/-  (1.09e-105, 4.67e-53j)
First non-zero Taylor coeff = 0.79212283864603056936 +/- 6.0749e-52
L(6.5) = (0.83934551203194208649 - 2.1002674689610803795e-51j)  +/-  (1.29e-37, 1.29e-37j)
L(7.5) = (0.90737569627003168218 - 1.5111677113148223163e-41j)  +/-  (3.4e-36, 3.4e-36j)
First 10 zeros
Zero 0 = 9.2223793999211025222 +/- 2.1895e-47
Zero 1 = 13.907549861392134406 +/- 7.0065e-46
Zero 2 = 17.442776978234473314 +/- 8.9683e-44
Zero 3 = 19.656513141954961 +/- 3.5873e-43
Zero 4 = 22.336103637209867276 +/- 1.8367e-40
Zero 5 = 25.274636548112365357 +/- 9.404e-38
Zero 6 = 26.804391158350403033 +/- 1.1755e-38
Zero 7 = 28.831682624186875445 +/- 9.404e-38
Zero 8 = 31.178209498360259064 +/- 6.0185e-36
Zero 9 = 32.774875382231207442 +/- 1.9259e-34
Z-plot in [0, 10]:
0.00	0.79	                              |----o
0.50	0.80	                              |----o
1.00	0.82	                              |-----o
1.50	0.85	                              |-----o
2.00	0.90	                              |-----o
2.50	0.96	                              |------o
3.00	1.03	                              |------o
3.50	1.12	                              |-------o
4.00	1.20	                              |--------o
4.50	1.29	                              |--------o
5.00	1.37	                              |---------o
5.50	1.44	                              |---------o
6.00	1.47	                              |----------o
6.50	1.46	                              |---------o
7.00	1.39	                              |---------o
7.50	1.24	                              |--------o
8.00	1.00	                              |------o
8.50	0.66	                              |---o
9.00	0.22	                              |o
9.22	zero	                              Z
9.50	-0.29	                            o-|
10.00	-0.84	                        o-----|
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
// https://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/1/12/a/a/
map<int64_t, vector<int64_t>> euler_factors  = {
 {2, {1, 24, 2048}},
 {3, {1, -252, 177147}},
 {5, {1, -4830, 48828125}},
 {7, {1, 16744, 1977326743}},
 {11, {1, -534612, 285311670611}},
 {13, {1, 577738, 1792160394037}},
 {17, {1, 6905934, 34271896307633}},
 {19, {1, -10661420, 116490258898219}},
 {23, {1, -18643272, 952809757913927}},
 {29, {1, -128406630, 12200509765705829}},
 {31, {1, 52843168, 25408476896404831}},
 {37, {1, 182213314, 177917621779460413}},
 {41, {1, -308120442, 550329031716248441}},
 {43, {1, 17125708, 929293739471222707}},
 {47, {1, -2687348496, 2472159215084012303}},
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

  // we have a degree 2 L-function with motivic weight 11
  L = Lfunc_init(2, 1, 5.5, mus, &ecode);
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
  printf("Epsilon = ");
  acb_printd(Lfunc_epsilon(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_epsilon(L) << endl;
  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_Taylor(L) << endl;


  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 6.5, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(6.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  ecode|=Lfunc_special_value(ctmp, L, 7.5,0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  printf("L(7.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
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
