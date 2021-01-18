// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Make up a degree 2 L-function associated the Elliptic Curve with LMFDB label 389.a1
 * https://www.lmfdb.org/EllipticCurve/Q/389/a/1
 *
 * Change the following two lines, to modify the number of decimal digits printed, or to print raw format (not human friendly)
 */
#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file should generate something
 * comparable to:
Order of vanishing = 2
Epsilon = (1 + 0j)  +/-  (4.13e-115, 9.09e-58j)
First non-zero Taylor coeff = 0.75931650028842677023 +/- 2.115e-51
L(1.5) = (0.13376843254631841442 - 1.7093870164286589522e-57j)  +/-  (7.9e-40, 7.9e-40j)
L(2.5) = (0.55297586704645019242 - 2.2906285188741987293e-48j)  +/-  (9.35e-42, 9.35e-42j)
First 10 zeros
Zero 0 = 2.8760990712604652018 +/- 1.044e-53
Zero 1 = 4.4168960836652578292 +/- 2.0881e-53
Zero 2 = 5.7934026339283652715 +/- 2.0881e-53
Zero 3 = 6.985966652828689218 +/- 6.6819e-52
Zero 4 = 7.4749074957854308894 +/- 6.6819e-52
Zero 5 = 8.6332052445633262416 +/- 2.6728e-51
Zero 6 = 9.6330788021849134547 +/- 2.1382e-50
Zero 7 = 10.351433312881496671 +/- 3.4211e-49
Zero 8 = 11.110935538806796228 +/- 3.4211e-49
Zero 9 = 11.933527327884206155 +/- 1.7106e-49
Z-plot in [0, 10]:
0.00	-0.00	                              o
0.50	0.21	                              |o
1.00	0.96	                              |------o
1.50	2.22	                              |---------------o
2.00	3.12	                              |----------------------o
2.50	2.21	                              |---------------o
2.88	zero	                              Z
3.00	-0.90	                        o-----|
3.50	-3.94	 o----------------------------|
4.00	-3.44	     o------------------------|
4.42	zero	                              Z
4.50	0.80	                              |-----o
5.00	4.07	                              |-----------------------------o
5.50	2.40	                              |-----------------o
5.79	zero	                              Z
6.00	-1.44	                    o---------|
6.50	-2.01	               o--------------|
6.99	zero	                              Z
7.00	0.05	                              o
7.47	zero	                              Z
7.50	-0.08	                              o
8.00	-1.79	                 o------------|
8.50	-0.71	                         o----|
8.63	zero	                              Z
9.00	1.50	                              |----------o
9.50	0.57	                              |---o
9.63	zero	                              Z
10.00	-0.80	                        o-----|
10.35	zero	                              Z
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
  {3, {1, 2, 3}},
  {5, {1, 3, 5}},
  {7, {1, 5, 7}},
  {11, {1, 4, 11}},
  {13, {1, 3, 13}},
  {17, {1, 6, 17}},
  {19, {1, -5, 19}},
  {23, {1, 4, 23}},
  {29, {1, 6, 29}},
  {31, {1, -4, 31}},
  {37, {1, 8, 37}},
  {41, {1, 3, 41}},
  {43, {1, -12, 43}},
  {47, {1, 2, 47}},
  {53, {1, 6, 53}},
  {59, {1, -3, 59}},
  {61, {1, 8, 61}},
  {67, {1, 5, 67}},
  {71, {1, 10, 71}},
  {73, {1, 7, 73}},
  {79, {1, 13, 79}},
  {83, {1, 12, 83}},
  {89, {1, 8, 89}},
  {97, {1, 9, 97}},
  {101, {1, 4, 101}},
  {103, {1, 6, 103}},
  {107, {1, 0, 107}},
  {109, {1, -2, 109}},
  {113, {1, -9, 113}},
  {127, {1, -11, 127}},
  {131, {1, -2, 131}},
  {137, {1, 2, 137}},
  {139, {1, -16, 139}},
  {149, {1, -10, 149}},
  {151, {1, 10, 151}},
  {157, {1, 6, 157}},
  {163, {1, 20, 163}},
  {167, {1, -8, 167}},
  {173, {1, -15, 173}},
  {179, {1, 19, 179}},
  {181, {1, -15, 181}},
  {191, {1, 4, 191}},
  {193, {1, 18, 193}},
  {197, {1, -8, 197}},
  {199, {1, 2, 199}},
  {211, {1, -7, 211}},
  {223, {1, 5, 223}},
  {227, {1, -12, 227}},
  {229, {1, 8, 229}},
  {233, {1, 6, 233}},
  {239, {1, -5, 239}},
  {241, {1, 30, 241}},
  {251, {1, -22, 251}},
  {257, {1, -8, 257}},
  {263, {1, 24, 263}},
  {269, {1, -3, 269}},
  {271, {1, 8, 271}},
  {277, {1, -17, 277}},
  {281, {1, 12, 281}},
  {283, {1, 9, 283}},
  {293, {1, 10, 293}},
  {307, {1, -26, 307}},
  {311, {1, 7, 311}},
  {313, {1, -19, 313}},
  {317, {1, 22, 317}},
  {331, {1, -4, 331}},
  {337, {1, 18, 337}},
  {347, {1, -9, 347}},
  {349, {1, 10, 349}},
  {353, {1, -14, 353}},
  {359, {1, 20, 359}},
  {367, {1, -2, 367}},
  {373, {1, 22, 373}},
  {379, {1, -4, 379}},
  {383, {1, 7, 383}},
  {389, {1, -1}},
  {397, {1, -24, 397}},
  {401, {1, 22, 401}},
  {409, {1, 1, 409}},
  {419, {1, 35, 419}},
  {421, {1, 16, 421}},
  {431, {1, -11, 431}},
  {433, {1, -29, 433}},
  {439, {1, -16, 439}},
  {443, {1, -35, 443}},
  {449, {1, 24, 449}},
  {457, {1, 9, 457}},
  {461, {1, 14, 461}}
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
  L = Lfunc_init(2, 389, 0.5, mus, &ecode);
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
  char bsd_str[] = "0.759316500288426770230192607894722019078097516494924351585805092547991072747780552440215935002328859166877149206819857682526013110483564490284990627405";
  arb_init(bsd);
  arb_set_str(bsd, bsd_str, 400);
  assert(arb_overlaps(Lfunc_Taylor(L), bsd));
  arb_clear(bsd);

  //FIXME add asserts and these two lines to intro text
  // ~ 0.133768432546318414418290373661
  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 1.5, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  // ~ 0.552975867046450192416260240311
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
