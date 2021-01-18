// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Make up a degree 2 L-function associated to the classical modular form 23.1.b.a
 * http://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/23/1/b/a/
 *
 * Change the following two lines, to modify the number of decimal digits printed, or to print raw format (not human friendly)
 */
#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file should generate something
 * comparable to:
Order of vanishing = 0
Epsilon = (1 + 0j)  +/-  (1.1e-117, 4.7e-59j)
First non-zero Taylor coeff = 0.1740363269879341835 +/- 8.2317e-59
L(1) = (0.36840932071582682111 - 1.2462241907180155708e-58j)  +/-  (3.26e-39, 3.26e-39j)
L(2) = (0.6747996946478415583 - 3.4460094331677626312e-48j)  +/-  (1.44e-42, 1.44e-42j)
First 10 zeros
Zero 0 = 5.1156833288151175986 +/- 4.1762e-53
Zero 1 = 7.1592622905417038413 +/- 1.6705e-52
Zero 2 = 8.881396573689177475 +/- 2.6728e-51
Zero 3 = 10.282027400852132054 +/- 2.1382e-50
Zero 4 = 11.430036353048133119 +/- 4.2764e-50
Zero 5 = 12.934409667718411836 +/- 1.7106e-49
Zero 6 = 14.662480954574659126 +/- 5.4738e-48
Zero 7 = 16.498232513362401342 +/- 8.7581e-47
Zero 8 = 17.101263303581831552 +/- 7.0065e-46
Zero 9 = 18.080661743444214077 +/- 1.4013e-45
Z-plot in [0, 10]:
0.00	0.17	                              |o
0.50	0.26	                              |o
1.00	0.50	                              |--o
1.50	0.88	                              |-----o
2.00	1.41	                              |---------o
2.50	2.05	                              |--------------o
3.00	2.64	                              |------------------o
3.50	2.95	                              |---------------------o
4.00	2.73	                              |-------------------o
4.50	1.82	                              |------------o
5.00	0.37	                              |-o
5.12	zero	                              Z
5.50	-1.15	                      o-------|
6.00	-2.04	               o--------------|
6.50	-1.79	                 o------------|
7.00	-0.51	                           o--|
7.16	zero	                              Z
7.50	1.01	                              |------o
8.00	1.71	                              |-----------o
8.50	1.07	                              |-------o
8.88	zero	                              Z
9.00	-0.34	                            o-|
9.50	-1.21	                     o--------|
10.00	-0.72	                         o----|
10.28	zero	                              Z
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
  {2, {1, 1, 1}},
  {3, {1, 1, 1}},
  {5, {1, 0, -1}},
  {7, {1, 0, -1}},
  {11, {1, 0, -1}},
  {13, {1, 1, 1}},
  {17, {1, 0, -1}},
  {19, {1, 0, -1}},
  {23, {1, -1}},
  {29, {1, 1, 1}},
  {31, {1, 1, 1}},
  {37, {1, 0, -1}},
  {41, {1, 1, 1}},
  {43, {1, 0, -1}},
  {47, {1, 1, 1}},
  {53, {1, 0, -1}},
  {59, {1, -2, 1}},
  {61, {1, 0, -1}},
  {67, {1, 0, -1}},
  {71, {1, 1, 1}},
  {73, {1, 1, 1}},
  {79, {1, 0, -1}},
  {83, {1, 0, -1}},
  {89, {1, 0, -1}},
  {97, {1, 0, -1}},
  {101, {1, -2, 1}},
  {103, {1, 0, -1}},
  {107, {1, 0, -1}},
  {109, {1, 0, -1}},
  {113, {1, 0, -1}}
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

  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  L = Lfunc_init(2, 23, 0.0, mus, &ecode);
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


  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 1.0, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  ecode|=Lfunc_special_value(ctmp, L, 2,0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  printf("L(2) = ");acb_printd(ctmp, DIGITS);printf("\n");
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
