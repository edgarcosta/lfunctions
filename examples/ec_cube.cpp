// Copyright Edgar Costa 2024
// See LICENSE file for license details.
/*
 * Computes the 3rd power of the Elliptic Curve L-function 37.a1: L(E, s)^3
 * Demonstrates the power extraction feature (extract_powers = YES).
 *
 * It supplies the exact local factors L_p(E, T)^3. The library's power guard
 * detects that the 2nd moment corresponds to a 3rd power and that the conductor
 * is a perfect cube, securely extracts the underlying E, computes it at 
 * lower degree/conductor, and assembles the final result.
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
#include "glfunc.h"
#include "glfunc_internals.h"
#include <cassert>

using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;

#define DIGITS 20

// LMFDB label 37.a1
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

// Supply L_p(T)^3
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  acb_poly_zero(poly);
  auto it = euler_factors.find(p);
  if( it != euler_factors.end() ) {
    acb_poly_t ep;
    acb_poly_init(ep);
    for(size_t i = 0; i < it->second.size(); ++i)
      acb_poly_set_coeff_si(ep, i, it->second[i]);
    
    acb_poly_pow_ui(poly, ep, 3, prec);
    acb_poly_clear(ep);
  }
}

int main ()
{
  Lerror_t ecode = ERR_SUCCESS;

  double mus[6] = {0, 0, 0, 1, 1, 1};
  Lparams_t Lp;
  Lp.degree = 6;
  Lp.conductor = 37ull * 37ull * 37ull; // 50653
  Lp.normalisation = 0.5;
  Lp.mus = mus;
  Lp.target_prec = 100;
  Lp.wprec = 0;
  Lp.gprec = 0;
  Lp.self_dual = YES;
  Lp.rank = DK;
  Lp.cache_dir = (char *)".";
  Lp.extract_powers = YES; // enable power extraction

  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    return 0;
  }

  // populate local factors
  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    return 0;
  }

  // do the computation
  ecode |= Lfunc_compute(L);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    return 0;
  }

  printf("L-function successfully computed via extraction!\n");
  printf("Rank = %" PRIu64 " (which is 3 * 1)\n", Lfunc_rank(L));
  printf("Sign = ");
  acb_printd(Lfunc_sign(L), DIGITS);
  printf(" (which is (-1)^3)\n");

  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L), DIGITS);
  printf("\n");

  printf("First isolated zero (repeats 3 times):\n");
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  for(int i = 0; i < 4; ++i) {
    printf("Zero %d = ", i);
    arb_printd(zeros + i, DIGITS);
    printf("\n");
  }

  acb_t ctmp; acb_init(ctmp);
  ecode |= Lfunc_special_value(ctmp, L, 1.5, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  printf("L(1.5) = "); acb_printd(ctmp, DIGITS); printf("\n");
  acb_clear(ctmp);

  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr, ecode);

  return 0;
}
