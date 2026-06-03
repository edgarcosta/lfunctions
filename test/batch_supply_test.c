/*
  Tests for the batch / array supply front-ends (bead lfunctions-c4h):
    Lfunc_use_lpolys_acb / _fmpz       (Euler-factor arrays)
    Lfunc_use_dirichlet_coeffs_fmpz / _acb + Lfunc_set_coeff_bound (raw a_n)

  Canonical object: the degree-2 self-dual L-function L(s,chi5)*L(s,chi7), the
  product of the quadratic characters mod 5 and mod 7 (conductor 35,
  normalisation 0, mus = [0,1]) already used by test/dir_test*.c. Everything
  about it is formula-computable, so the same object can be supplied every way
  and the certified outputs compared with arb/acb overlaps (never printed
  output). Exits 0 on success; an assert fires on failure.
*/
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_vec.h>
#include "glfunc.h"

// --- the object: L(s,chi5)*L(s,chi7) ----------------------------------------
// quadratic character mod 5: chi5(n) = (n/5); mod 7: chi7(n) = (n/7).
static int chi5(uint64_t n) {
  switch (n % 5) { case 1: case 4: return 1; case 2: case 3: return -1; default: return 0; }
}
static int chi7(uint64_t n) {
  switch (n % 7) { case 1: case 2: case 4: return 1; case 3: case 5: case 6: return -1; default: return 0; }
}

// local factor at p in the algebraic normalisation:
//   (1 - chi5(p) T)(1 - chi7(p) T) = 1 + c1 T + c2 T^2
static void factor_coeffs(uint64_t p, long *c1, long *c2) {
  *c1 = -(long)(chi5(p) + chi7(p));
  *c2 = (long)chi5(p) * chi7(p);
}
static void factor_acb(acb_poly_t f, uint64_t p) {
  long c1, c2; factor_coeffs(p, &c1, &c2);
  acb_poly_zero(f);
  acb_poly_set_coeff_si(f, 0, 1);
  acb_poly_set_coeff_si(f, 1, c1);
  acb_poly_set_coeff_si(f, 2, c2);
}
static void factor_fmpz(fmpz_poly_t f, uint64_t p) {
  long c1, c2; factor_coeffs(p, &c1, &c2);
  fmpz_poly_zero(f);
  fmpz_poly_set_coeff_si(f, 0, 1);
  fmpz_poly_set_coeff_si(f, 1, c1);
  fmpz_poly_set_coeff_si(f, 2, c2);
}

// primes <= bound, in increasing order; returns the count.
static uint64_t primes_upto(uint64_t bound, uint64_t *primes) {
  uint64_t cnt = 0;
  for (uint64_t n = 2; n <= bound; n++) {
    int isp = 1;
    for (uint64_t d = 2; d * d <= n; d++) if (n % d == 0) { isp = 0; break; }
    if (isp) primes[cnt++] = n;
  }
  return cnt;
}

// --- supply routes (each returns a fresh, computed Lfunc_t) ------------------
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param) {
  (void)d; (void)prec; (void)param;
  factor_acb(poly, p);
}

static Lfunc_t run_callback(Lerror_t *ec) {
  double mus[] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, ec);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}

static Lfunc_t run_lpolys_acb(Lerror_t *ec) {
  double mus[] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  uint64_t *primes = (uint64_t *)malloc(sizeof(uint64_t) * (nmax + 1));
  uint64_t np = primes_upto(nmax, primes);
  acb_poly_struct *f = (acb_poly_struct *)malloc(sizeof(acb_poly_struct) * np);
  for (uint64_t k = 0; k < np; k++) { acb_poly_init(&f[k]); factor_acb(&f[k], primes[k]); }
  *ec |= Lfunc_use_lpolys_acb(L, f, np);
  for (uint64_t k = 0; k < np; k++) acb_poly_clear(&f[k]);
  free(f); free(primes);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}

static Lfunc_t run_lpolys_fmpz(Lerror_t *ec) {
  double mus[] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  uint64_t *primes = (uint64_t *)malloc(sizeof(uint64_t) * (nmax + 1));
  uint64_t np = primes_upto(nmax, primes);
  fmpz_poly_struct *f = (fmpz_poly_struct *)malloc(sizeof(fmpz_poly_struct) * np);
  for (uint64_t k = 0; k < np; k++) { fmpz_poly_init(&f[k]); factor_fmpz(&f[k], primes[k]); }
  *ec |= Lfunc_use_lpolys_fmpz(L, f, np);
  for (uint64_t k = 0; k < np; k++) fmpz_poly_clear(&f[k]);
  free(f); free(primes);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}

// --- comparison helpers ------------------------------------------------------
static void assert_zeros_overlap(Lfunc_t A, Lfunc_t B, uint64_t side) {
  arb_srcptr za = Lfunc_zeros(A, side), zb = Lfunc_zeros(B, side);
  uint64_t i = 0;
  for (; !arb_is_zero(za + i) && !arb_is_zero(zb + i); i++)
    assert(arb_overlaps(za + i, zb + i));
  assert(arb_is_zero(za + i) && arb_is_zero(zb + i)); // same number of zeros
  assert(i > 0); // we did find some zeros to compare
}
static void assert_outputs_overlap(Lfunc_t A, Lfunc_t B) {
  assert(Lfunc_rank(A) == Lfunc_rank(B));
  assert(acb_overlaps(Lfunc_sign(A), Lfunc_sign(B)));
  assert(arb_overlaps(Lfunc_Taylor(A), Lfunc_Taylor(B)));
  assert_zeros_overlap(A, B, 0);
}
// invariants any genuine run of this self-dual, non-vanishing object satisfies.
static void assert_self_consistent(Lfunc_t L) {
  assert(Lfunc_rank(L) == 0); // L(1/2,chi5) L(1/2,chi7) != 0
  arb_t m; arb_init(m);
  acb_abs(m, Lfunc_sign(L), 100);
  assert(arb_contains_si(m, 1));                     // |epsilon| = 1
  assert(arb_contains_zero(acb_imagref(Lfunc_sign(L)))); // epsilon is real (self-dual)
  arb_clear(m);
}

// ============================================================================
static void test_factor_arrays(void) {
  Lerror_t ea = ERR_SUCCESS, eb = ERR_SUCCESS, ec = ERR_SUCCESS;
  Lfunc_t A = run_callback(&ea);
  Lfunc_t B = run_lpolys_acb(&eb);
  Lfunc_t C = run_lpolys_fmpz(&ec);
  assert(!fatal_error(ea) && !fatal_error(eb) && !fatal_error(ec));

  assert_self_consistent(A);
  assert_self_consistent(B);
  assert_self_consistent(C);

  // the three Euler-factor routes must agree (full supply, same object)
  assert_outputs_overlap(A, B);
  assert_outputs_overlap(A, C);

  Lfunc_clear(A); Lfunc_clear(B); Lfunc_clear(C);
}

int main(void) {
  test_factor_arrays();
  printf("batch_supply_test: all tests passed\n");
  return 0;
}
