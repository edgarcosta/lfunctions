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
#include <cstdio>
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
#include <flint/fmpz_poly.h>
#include "glfunc.h"
#include <cassert>

using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;

#define DIGITS 20

static int report_fatal(Lerror_t ecode, Lfunc_t L)
{
  fprint_errors(stderr, ecode);
  if(L)
    Lfunc_clear(L);
  return 1;
}

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

// Supply exact L_p(T)^3
void lpoly_callback(fmpz_poly_t poly, uint64_t p, int d __attribute__((unused)), void *param __attribute__((unused)))
{
  fmpz_poly_zero(poly);
  auto it = euler_factors.find(p);
  if( it != euler_factors.end() ) {
    fmpz_poly_t ep;
    fmpz_poly_init(ep);
    for(size_t i = 0; i < it->second.size(); ++i)
      fmpz_poly_set_coeff_si(ep, i, it->second[i]);
    
    fmpz_poly_pow(poly, ep, 3);
    fmpz_poly_clear(ep);
  }
}

int main ()
{
  Lerror_t ecode = ERR_SUCCESS;

  double mus[6] = {0, 0, 0, 1, 1, 1};
  Lparams_t Lp = {};
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
    return report_fatal(ecode, NULL);
  }

  // populate local factors
  ecode |= Lfunc_use_all_lpolys_fmpz(L, lpoly_callback, NULL);
  if(fatal_error(ecode)) {
    return report_fatal(ecode, L);
  }

  // do the computation
  ecode |= Lfunc_compute(L);
  if(fatal_error(ecode)) {
    return report_fatal(ecode, L);
  }

  Lfunc_t *factors = NULL;
  uint64_t *mults = NULL;
  uint64_t nfactors = Lfunc_factors(L, &factors, &mults);
  assert(nfactors == 1);
  assert(mults[0] == 3);
  assert(factors[0] != NULL);
  Lfunc_t E = factors[0];

  int64_t factor_rank = Lfunc_rank(E);
  int64_t rank = Lfunc_rank(L);
  assert(factor_rank == 1);
  assert(rank == 3 * factor_rank);

  acb_t sign_cube;
  acb_init(sign_cube);
  acb_pow_ui(sign_cube, Lfunc_sign(E), 3, Lfunc_wprec(L));
  assert(acb_overlaps(sign_cube, Lfunc_sign(L)));
  acb_clear(sign_cube);

  arb_srcptr factor_zeros = Lfunc_zeros(E, 0);
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  assert(arb_overlaps(zeros + 0, factor_zeros + 0));
  assert(arb_overlaps(zeros + 1, factor_zeros + 0));
  assert(arb_overlaps(zeros + 2, factor_zeros + 0));
  assert(arb_overlaps(zeros + 3, factor_zeros + 1));

  arb_t taylor_cube;
  arb_init(taylor_cube);
  arb_pow_ui(taylor_cube, Lfunc_Taylor(E), 3, Lfunc_wprec(L));
  assert(arb_overlaps(Lfunc_Taylor(L), taylor_cube));
  arb_clear(taylor_cube);

  printf("L-function successfully computed via extraction!\n");
  printf("Extracted factors = %" PRIu64 ", multiplicity = %" PRIu64 "\n", nfactors, mults[0]);
  printf("Rank = %" PRId64 " (which is 3 * %" PRId64 ")\n", rank, factor_rank);
  printf("Sign = ");
  acb_printd(Lfunc_sign(L), DIGITS);
  printf(" (which is (-1)^3)\n");

  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L), DIGITS);
  printf("\n");

  printf("First isolated zero (repeats 3 times):\n");
  for(int i = 0; i < 4; ++i) {
    printf("Zero %d = ", i);
    arb_printd(zeros + i, DIGITS);
    printf("\n");
  }

  acb_t ctmp, etmp, ecube;
  acb_init(ctmp); acb_init(etmp); acb_init(ecube);
  Lerror_t sv = ERR_SUCCESS;
  sv |= Lfunc_special_value(etmp, E, 1.5, 0.0);
  sv |= Lfunc_special_value(ctmp, L, 1.5, 0.0);
  ecode |= sv;
  if(fatal_error(sv)) {
    acb_clear(ctmp); acb_clear(etmp); acb_clear(ecube);
    return report_fatal(sv, L);
  }
  acb_pow_ui(ecube, etmp, 3, Lfunc_wprec(L));
  assert(acb_overlaps(ctmp, ecube));
  printf("L(1.5) = "); acb_printd(ctmp, DIGITS); printf("\n");
  acb_clear(ctmp); acb_clear(etmp); acb_clear(ecube);

  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr, ecode);

  return 0;
}
