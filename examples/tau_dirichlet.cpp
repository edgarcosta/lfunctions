// Copyright Edgar Costa 2026
// See LICENSE file for license details.
/*
 * The Ramanujan tau L-function (degree 2, conductor 1, motivic weight 11),
 * supplied via the raw Dirichlet-coefficient front-end instead of Euler factors.
 *
 * examples/tau.cpp hands over the local factors 1 - tau(p) T + p^11 T^2 and lets
 * the library build the Dirichlet coefficients. Here we instead hand over the
 * coefficients tau(n) directly with Lfunc_use_dirichlet_coeffs_fmpz, which is the
 * natural entry point when you already hold a_n (say from a database row). The
 * library validates each tau(n) against the degree-2 Euler-product bound
 * automatically (no caller-supplied bound needed). The default Turing verifier
 * can still run from raw coefficients; Buthe-only builds lack the per-prime data
 * for RH verification. The table currently reaches Lfunc_nmax; if a future
 * precision/bound change outgrows it, the example will surface ERR_INSUFF_EULER.
 *
 * See https://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/1/12/a/a/.
 */
#define __STDC_FORMAT_MACROS
#include <cassert>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <flint/acb.h>
#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include "glfunc.h"

using std::int64_t;
using std::map;
using std::size_t;
using std::vector;

// Local factors 1 - tau(p) T + p^11 T^2, the same data tau.cpp uses (LMFDB).
static const map<int64_t, vector<int64_t>> euler_factors = {
  {2,  {1, 24, 2048}},
  {3,  {1, -252, 177147}},
  {5,  {1, -4830, 48828125}},
  {7,  {1, 16744, 1977326743}},
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

// tau(p^e) from tau(p) = -c_1 and p^11 = c_2 via the Hecke recurrence
//   tau(p^k) = tau(p) tau(p^{k-1}) - p^11 tau(p^{k-2}).
static int64_t tau_prime_power(int64_t p, int e) {
  const vector<int64_t> &f = euler_factors.at(p);
  int64_t taup = -f[1], p11 = f[2];
  int64_t a = 1, b = taup; // tau(p^0), tau(p^1)
  if (e == 0) return a;
  for (int k = 2; k <= e; ++k) { int64_t c = taup * b - p11 * a; a = b; b = c; }
  return b;
}

// tau(n) for n = 1..N by multiplicativity. Needs tau(p) for every prime <= N,
// so N must stay within the table above (primes up to 47, i.e. n <= 52).
static void build_tau(vector<int64_t> &tau, size_t N) {
  tau.assign(N + 1, 0);
  if (N >= 1) tau[1] = 1;
  for (size_t n = 2; n <= N; ++n) {
    int64_t p = 0;
    for (const auto &kv : euler_factors)
      if (n % (size_t)kv.first == 0) { p = kv.first; break; }
    assert(p && "tau table does not cover a prime factor of n");
    int e = 0; size_t m = n;
    while (m % (size_t)p == 0) { m /= (size_t)p; ++e; }
    tau[n] = tau_prime_power(p, e) * tau[m];
  }
}

int main() {
  Lfunc_t L;
  double mus[2] = {0, 1};
  Lerror_t ecode;

  // degree 2, conductor 1, motivic weight 11 -> normalisation 11/2 = 5.5
  L = Lfunc_init(2, 1, 5.5, mus, &ecode);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

  // tau(p) is tabulated up to p = 47, so we can build tau(n) for n <= 52. This
  // currently reaches nmax; if the library wants more after a future bound
  // change, raw coefficient supply reduces nmax and warns.
  const size_t TAU_MAX = 52;
  uint64_t nmax = Lfunc_nmax(L);
  size_t len = nmax < TAU_MAX ? (size_t)nmax : TAU_MAX;
  vector<int64_t> tau;
  build_tau(tau, len);

  // Hand over the coefficients directly. LFUNC_ALGEBRAIC_NORM: these are the algebraic
  // tau(n); the library applies the n^{-5.5} shift to the analytic normalisation
  // and validates |tau(n)| n^{-11/2} <= d(n) <= n against the degree-2 bound.
  fmpz *a = _fmpz_vec_init(len);
  for (size_t n = 1; n <= len; ++n) fmpz_set_si(a + (n - 1), tau[n]);
  ecode |= Lfunc_use_dirichlet_coeffs_fmpz(L, a, len, LFUNC_ALGEBRAIC_NORM);
  _fmpz_vec_clear(a, len);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); Lfunc_clear(L); return 1; }

  ecode |= Lfunc_compute(L);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); Lfunc_clear(L); return 1; }

  // In the default TURING build, raw coefficients still run the RH check;
  // Buthe-only builds report unavailable. The short-supply warning appears only
  // if the current nmax outgrows the table above.
  if (len < nmax)
    assert(ecode & ERR_INSUFF_EULER);
  else
    assert(!(ecode & ERR_INSUFF_EULER));
#ifdef TURING
  assert(!(ecode & ERR_RH_UNAVAILABLE));
#else
  assert(ecode & ERR_RH_UNAVAILABLE);
#endif
  assert(!fatal_error(ecode));

  printf("Order of vanishing = %" PRId64 "\n", Lfunc_rank(L));
  printf("Sign = "); acb_printd(Lfunc_sign(L), 20); printf("\n");
  printf("First non-zero Taylor coeff = "); arb_printd(Lfunc_Taylor(L), 20); printf("\n");

  acb_t ctmp; acb_init(ctmp);
  ecode |= Lfunc_special_value(ctmp, L, 6.5, 0.0);
  if (fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    acb_clear(ctmp);
    Lfunc_clear(L);
    return 1;
  }
  printf("L(6.5) = "); acb_printd(ctmp, 20); printf("\n");
  { // same computed value as tau.cpp's factor route
    acb_t ref; acb_init(ref);
    arb_set_str(acb_realref(ref), "0.83934551203194208649", 300);
    arb_set_str(acb_imagref(ref), "0", 300);
    arb_add_error_2exp_si(acb_realref(ref), -50);
    arb_add_error_2exp_si(acb_imagref(ref), -50);
    assert(acb_overlaps(ctmp, ref));
    acb_clear(ref);
  }
  acb_clear(ctmp);

  arb_srcptr zeros = Lfunc_zeros(L, 0);
  printf("First zero = "); arb_printd(zeros + 0, 20); printf("\n");
  { // same computed first zero ordinate as tau.cpp
    arb_t ref; arb_init(ref);
    arb_set_str(ref, "9.2223793999211025222", 300);
    arb_add_error_2exp_si(ref, -50);
    assert(arb_overlaps(zeros + 0, ref));
    arb_clear(ref);
  }

  Lfunc_clear(L);
  fprint_errors(stderr, ecode); // prints the shortfall/RH warning notices
  return 0;
}
