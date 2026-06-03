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

// raw Dirichlet coefficient a_n = sum_{d|n} chi5(d) chi7(n/d) (the convolution
// whose Dirichlet series is exactly prod_p (1-chi5(p)T)(1-chi7(p)T) inverted).
static long an(uint64_t n) {
  long s = 0;
  for (uint64_t d = 1; d <= n; d++)
    if (n % d == 0) s += (long)chi5(d) * chi7(n / d);
  return s;
}

// the degree-2 default bound coeff_bound(2) = 1, alpha = 1: |a_n| <= d(n) <= n.
static void set_default_bound(Lfunc_t L, Lerror_t *ec) {
  arb_t C; arb_init(C); arb_set_ui(C, 1);
  *ec |= Lfunc_set_coeff_bound(L, C, 1.0);
  arb_clear(C);
}

static Lfunc_t run_raw_fmpz(Lerror_t *ec) {
  double mus[] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  set_default_bound(L, ec);
  fmpz *a = _fmpz_vec_init(nmax);
  for (uint64_t n = 1; n <= nmax; n++) fmpz_set_si(a + (n - 1), an(n));
  *ec |= Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM); // norm 0: no shift
  _fmpz_vec_clear(a, nmax);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}

static Lfunc_t run_raw_acb(Lerror_t *ec) {
  double mus[] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  set_default_bound(L, ec);
  acb_ptr a = _acb_vec_init(nmax);
  for (uint64_t n = 1; n <= nmax; n++) acb_set_si(a + (n - 1), an(n));
  *ec |= Lfunc_use_dirichlet_coeffs_acb(L, a, nmax, ALGEBRAIC_NORM);
  _acb_vec_clear(a, nmax);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}

// Same analytic object as run_callback (conductor 35, analytic mus [0,1],
// analytic a_n = b_n) but parametrised with normalisation 0.5, mus [-0.5,0.5]
// (the [0,1],0.5 == [0.5,1.5],0 invariant). Supplying the algebraic numbers
// b_n * n^{0.5} under ALGEBRAIC_NORM must reproduce analytic b_n, agreeing with
// supplying b_n directly under ANALYTIC_NORM.
static Lfunc_t run_shift_algebraic(Lerror_t *ec) {
  double mus[] = {-0.5, 0.5};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.5, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  set_default_bound(L, ec);
  acb_ptr a = _acb_vec_init(nmax);
  arb_t rt; arb_init(rt);
  for (uint64_t n = 1; n <= nmax; n++) {
    arb_sqrt_ui(rt, n, 256);          // n^{0.5}
    arb_mul_si(rt, rt, an(n), 256);   // b_n * n^{0.5} (the algebraic coeff)
    acb_set_arb(a + (n - 1), rt);
  }
  arb_clear(rt);
  *ec |= Lfunc_use_dirichlet_coeffs_acb(L, a, nmax, ALGEBRAIC_NORM); // shift by n^{-0.5}
  _acb_vec_clear(a, nmax);
  if (fatal_error(*ec)) return L;
  *ec |= Lfunc_compute(L);
  return L;
}
static Lfunc_t run_shift_analytic(Lerror_t *ec) {
  double mus[] = {-0.5, 0.5};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.5, mus, ec);
  if (fatal_error(*ec)) return L;
  uint64_t nmax = Lfunc_nmax(L);
  set_default_bound(L, ec);
  acb_ptr a = _acb_vec_init(nmax);
  for (uint64_t n = 1; n <= nmax; n++) acb_set_si(a + (n - 1), an(n)); // analytic b_n
  *ec |= Lfunc_use_dirichlet_coeffs_acb(L, a, nmax, ANALYTIC_NORM);    // no shift
  _acb_vec_clear(a, nmax);
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

// raw a_n (fmpz) reproduces the callback's certified data.
static void test_raw_coeffs(void) {
  Lerror_t ea = ERR_SUCCESS, er = ERR_SUCCESS;
  Lfunc_t A = run_callback(&ea);
  Lfunc_t R = run_raw_fmpz(&er);
  assert(!fatal_error(ea) && !fatal_error(er));
  assert_self_consistent(R);
  assert_outputs_overlap(A, R);
  Lfunc_clear(A); Lfunc_clear(R);
}

// the fmpz and acb coefficient forms carry the same integers, so they must
// agree (both exact).
static void test_acb_vs_fmpz_coeffs(void) {
  Lerror_t ef = ERR_SUCCESS, eg = ERR_SUCCESS;
  Lfunc_t F = run_raw_fmpz(&ef);
  Lfunc_t G = run_raw_acb(&eg);
  assert(!fatal_error(ef) && !fatal_error(eg));
  assert_outputs_overlap(F, G);
  Lfunc_clear(F); Lfunc_clear(G);
}

// ALGEBRAIC_NORM with a nonzero normalisation (the library applies n^{-0.5})
// equals the same coefficients pre-shifted and supplied ANALYTIC_NORM, and both
// equal the canonical normalisation-0 run.
static void test_normalisation_flag(void) {
  Lerror_t ep = ERR_SUCCESS, eq = ERR_SUCCESS, ea = ERR_SUCCESS;
  Lfunc_t P = run_shift_algebraic(&ep);
  Lfunc_t Q = run_shift_analytic(&eq);
  Lfunc_t A = run_callback(&ea);
  assert(!fatal_error(ep) && !fatal_error(eq) && !fatal_error(ea));
  assert_self_consistent(P);
  assert_self_consistent(Q);
  assert_outputs_overlap(P, Q);
  assert_outputs_overlap(P, A);
  Lfunc_clear(P); Lfunc_clear(Q); Lfunc_clear(A);
}

// --- edge guards (fatal supply errors) --------------------------------------
static fmpz *make_an_fmpz(uint64_t len) {
  fmpz *a = _fmpz_vec_init(len);
  for (uint64_t n = 1; n <= len; n++) fmpz_set_si(a + (n - 1), an(n));
  return a;
}

// a_1 != 1 is rejected (fmpz: exact inequality).
static void test_guard_a1_fmpz(void) {
  double mus[] = {0, 1}; Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
  set_default_bound(L, &ec); assert(!fatal_error(ec));
  uint64_t nmax = Lfunc_nmax(L);
  fmpz *a = make_an_fmpz(nmax);
  fmpz_set_si(a + 0, 2); // a_1 = 2
  Lerror_t e = Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
  assert((e & ERR_A1_NOT_ONE) && fatal_error(e));
  _fmpz_vec_clear(a, nmax);
  Lfunc_clear(L);
}

// a_1 ball not containing 1 is rejected (acb).
static void test_guard_a1_acb(void) {
  double mus[] = {0, 1}; Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
  set_default_bound(L, &ec); assert(!fatal_error(ec));
  uint64_t nmax = Lfunc_nmax(L);
  acb_ptr a = _acb_vec_init(nmax);
  for (uint64_t n = 1; n <= nmax; n++) acb_set_si(a + (n - 1), an(n));
  acb_set_si(a + 0, 2); // a_1 = 2 exactly: does not contain 1
  Lerror_t e = Lfunc_use_dirichlet_coeffs_acb(L, a, nmax, ALGEBRAIC_NORM);
  assert((e & ERR_A1_NOT_ONE) && fatal_error(e));
  _acb_vec_clear(a, nmax);
  Lfunc_clear(L);
}

// raw a_n with no prior Lfunc_set_coeff_bound is rejected.
static void test_guard_missing_bound(void) {
  double mus[] = {0, 1}; Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
  uint64_t nmax = Lfunc_nmax(L);
  fmpz *a = make_an_fmpz(nmax); // a_1 = 1, valid; only the bound is missing
  Lerror_t e = Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
  assert((e & ERR_COEFF_BOUND) && fatal_error(e));
  _fmpz_vec_clear(a, nmax);
  Lfunc_clear(L);
}

// raw a_n incompatible with factors / a second raw supply.
static void test_guard_conflicts(void) {
  double mus[] = {0, 1};
  // (a) factor array, then raw -> conflict
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
    set_default_bound(L, &ec);
    uint64_t nmax = Lfunc_nmax(L);
    uint64_t *primes = (uint64_t *)malloc(sizeof(uint64_t) * (nmax + 1));
    uint64_t np = primes_upto(nmax, primes);
    acb_poly_struct *f = (acb_poly_struct *)malloc(sizeof(acb_poly_struct) * np);
    for (uint64_t k = 0; k < np; k++) { acb_poly_init(&f[k]); factor_acb(&f[k], primes[k]); }
    ec |= Lfunc_use_lpolys_acb(L, f, np);
    assert(!fatal_error(ec));
    for (uint64_t k = 0; k < np; k++) acb_poly_clear(&f[k]);
    free(f); free(primes);
    fmpz *a = make_an_fmpz(nmax);
    Lerror_t e = Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
    assert((e & ERR_SUPPLY_CONFLICT) && fatal_error(e));
    _fmpz_vec_clear(a, nmax);
    Lfunc_clear(L);
  }
  // (b) raw, then factor array -> conflict (either order)
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
    set_default_bound(L, &ec);
    uint64_t nmax = Lfunc_nmax(L);
    fmpz *a = make_an_fmpz(nmax);
    ec |= Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
    assert(!fatal_error(ec));
    _fmpz_vec_clear(a, nmax);
    uint64_t *primes = (uint64_t *)malloc(sizeof(uint64_t) * (nmax + 1));
    uint64_t np = primes_upto(nmax, primes);
    acb_poly_struct *f = (acb_poly_struct *)malloc(sizeof(acb_poly_struct) * np);
    for (uint64_t k = 0; k < np; k++) { acb_poly_init(&f[k]); factor_acb(&f[k], primes[k]); }
    Lerror_t e = Lfunc_use_lpolys_acb(L, f, np);
    assert((e & ERR_SUPPLY_CONFLICT) && fatal_error(e));
    for (uint64_t k = 0; k < np; k++) acb_poly_clear(&f[k]);
    free(f); free(primes);
    Lfunc_clear(L);
  }
  // (c) raw twice -> conflict
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
    set_default_bound(L, &ec);
    uint64_t nmax = Lfunc_nmax(L);
    fmpz *a = make_an_fmpz(nmax);
    ec |= Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
    assert(!fatal_error(ec));
    Lerror_t e = Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
    assert((e & ERR_SUPPLY_CONFLICT) && fatal_error(e));
    _fmpz_vec_clear(a, nmax);
    Lfunc_clear(L);
  }
}

// NEGATIVE soundness test: a declared bound the data violates must be caught,
// not silently accepted (otherwise M_error's tail bound is invalid). Declaring
// |a_n| <= 1 * n^0 = 1 is violated by a_9 = 3.
static void test_guard_violated_bound(void) {
  double mus[] = {0, 1}; Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec); assert(!fatal_error(ec));
  arb_t C; arb_init(C); arb_set_ui(C, 1);
  ec |= Lfunc_set_coeff_bound(L, C, 0.0); // alpha = 0 => bound is the constant 1
  arb_clear(C);
  assert(!fatal_error(ec));
  uint64_t nmax = Lfunc_nmax(L);
  assert(nmax >= 9); // a_9 = 3 must be in range to trip the check
  fmpz *a = make_an_fmpz(nmax); // a_1 = 1 (valid), but a_9 = 3 > 1
  Lerror_t e = Lfunc_use_dirichlet_coeffs_fmpz(L, a, nmax, ALGEBRAIC_NORM);
  assert((e & ERR_COEFF_BOUND) && fatal_error(e));
  _fmpz_vec_clear(a, nmax);
  Lfunc_clear(L);
}

int main(void) {
  test_factor_arrays();
  test_raw_coeffs();
  test_acb_vs_fmpz_coeffs();
  test_normalisation_flag();
  test_guard_a1_fmpz();
  test_guard_a1_acb();
  test_guard_missing_bound();
  test_guard_conflicts();
  test_guard_violated_bound();
  printf("batch_supply_test: all tests passed\n");
  return 0;
}
