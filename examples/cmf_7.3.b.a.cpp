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
Epsilon = (1 + 0j)  +/-  (2.52e-114, 2.24e-57j)
First non-zero Taylor coeff = 0.33298177148293397364 +/- 8.2499e-57
L(1) = (0.19673558706889585505 + 5.1688917606710579456e-57j)  +/-  (1.16e-38, 1.16e-38j)
L(2) = (0.46721176888437312394 - 1.2275191787387591289e-56j)  +/-  (2.75e-38, 2.75e-38j)
First 10 zeros
Zero 0 = 7.2145891812871844435 +/- 1.3364e-51
Zero 1 = 9.2568107485814985876 +/- 1.0691e-50
Zero 2 = 10.580974854111851196 +/- 1.3685e-48
Zero 3 = 12.965859683008166209 +/- 1.3685e-48
Zero 4 = 15.648104672594291642 +/- 7.0065e-46
Zero 5 = 16.796003155249041724 +/- 7.0065e-46
Zero 6 = 18.404976306799735748 +/- 1.121e-44
Zero 7 = 19.213914502107617614 +/- 8.9683e-44
Zero 8 = 20.72321683128359949 +/- 1.7937e-43
Zero 9 = 22.543649261523556672 +/- 1.1479e-41
Z-plot in [0, 10]:
0.00	0.33	                              |-o
0.50	0.36	                              |-o
1.00	0.45	                              |--o
1.50	0.61	                              |---o
2.00	0.85	                              |-----o
2.50	1.18	                              |-------o
3.00	1.58	                              |----------o
3.50	2.04	                              |--------------o
4.00	2.49	                              |-----------------o
4.50	2.86	                              |--------------------o
5.00	3.01	                              |---------------------o
5.50	2.85	                              |--------------------o
6.00	2.32	                              |----------------o
6.50	1.45	                              |---------o
7.00	0.43	                              |--o
7.21	zero	                              Z
7.50	-0.49	                           o--|
8.00	-1.00	                       o------|
8.50	-0.94	                       o------|
9.00	-0.37	                            o-|
9.26	zero	                              Z
9.50	0.32	                              |-o
10.00	0.63	                              |---o
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
{2, { 1, 3, 4 }},
{3, { 1, 0, -9 }},
{5, { 1, 0, -25 }},
{7, { 1, 7 }},
{11, { 1, 6, 121 }},
{13, { 1, 0, -169 }},
{17, { 1, 0, -289 }},
{19, { 1, 0, -361 }},
{23, { 1, -18, 529 }},
{29, { 1, 54, 841 }},
{31, { 1, 0, -961 }},
{37, { 1, 38, 1369 }},
{41, { 1, 0, -1681 }},
{43, { 1, -58, 1849 }},
{47, { 1, 0, -2209 }},
{53, { 1, 6, 2809 }},
{59, { 1, 0, -3481 }},
{61, { 1, 0, -3721 }},
{67, { 1, 118, 4489 }},
{71, { 1, -114, 5041 }},
{73, { 1, 0, -5329 }},
{79, { 1, 94, 6241 }},
{83, { 1, 0, -6889 }},
{89, { 1, 0, -7921 }},
{97, { 1, 0, -9409 }},
{101, { 1, 0, -10201 }},
{103, { 1, 0, -10609 }},
{107, { 1, -186, 11449 }},
{109, { 1, -106, 11881 }},
{113, { 1, 222, 12769 }},
{127, { 1, -2, 16129 }},
{131, { 1, 0, -17161 }},
{137, { 1, 174, 18769 }},
{139, { 1, 0, -19321 }},
{149, { 1, -186, 22201 }},
{151, { 1, -274, 22801 }},
{157, { 1, 0, -24649 }},
{163, { 1, -74, 26569 }},
{167, { 1, 0, -27889 }},
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

  // we have a degree 2 L-function with motivic weight 2
  L = Lfunc_init(2, 7, 1, mus, &ecode);
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
