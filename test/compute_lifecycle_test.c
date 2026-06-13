/*
  Regression test: Lfunc_compute lifecycle guards.

  Lfunc_compute divides the Dirichlet coefficients L->ans by sqrt(n) in place,
  so it is single-use, and it needs at least one supply route to have run (else
  M / M0 / dc were never set by Lfunc_nmax and the pipeline reads uninitialised
  state). Three misuse paths used to be unguarded:

    1. compute with no supply  -> read of uninitialised M/M0/dc (a segfault on
       the base);
    2. a second Lfunc_compute  -> re-divides ans by sqrt(n), silently changing a
       certified value (no fatal flag on the base);
    3. a supply call after compute -> pushes/overwrites already-normalised ans.

  All three are now fatal ERR_LIFECYCLE. Test object: the formula-computable
  degree-2 self-dual L(s,chi5)*L(s,chi7) (conductor 35, mus [0,1]) used by the
  other tests. Through the public glfunc.h API; exits 0 on success.
*/
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include "glfunc.h"

static int chi5(uint64_t n) {
  switch (n % 5) { case 1: case 4: return 1; case 2: case 3: return -1; default: return 0; }
}
static int chi7(uint64_t n) {
  switch (n % 7) { case 1: case 2: case 4: return 1; case 3: case 5: case 6: return -1; default: return 0; }
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param) {
  (void)d; (void)prec; (void)param;
  int c5 = chi5(p), c7 = chi7(p);
  acb_poly_zero(poly);
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_poly_set_coeff_si(poly, 1, -(c5 + c7));
  acb_poly_set_coeff_si(poly, 2, c5 * c7);
}

// (1) compute with no supply is fatal (and on the base it segfaulted reading
// uninitialised M/M0/dc).
static void test_compute_no_supply(void) {
  double mus[] = {0, 1};
  Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
  assert(!fatal_error(ec));
  Lerror_t e = Lfunc_compute(L);
  assert((e & ERR_LIFECYCLE) && fatal_error(e));
  Lfunc_clear(L);
}

// (2) a second compute is fatal AND the first certified result is left intact.
static void test_double_compute(void) {
  double mus[] = {0, 1};
  Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
  assert(!fatal_error(ec));
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  assert(!fatal_error(ec));
  ec |= Lfunc_compute(L);
  assert(!fatal_error(ec));

  arb_t taylor1; arb_init(taylor1);
  arb_set(taylor1, Lfunc_Taylor(L));   // certified leading coeff after compute 1
  assert(!arb_contains_zero(taylor1)); // rank 0: a genuine nonzero value

  Lerror_t e2 = Lfunc_compute(L);
  assert((e2 & ERR_LIFECYCLE) && fatal_error(e2)); // second compute rejected
  // The stored value must be untouched: on the base the re-division moved it.
  assert(arb_equal(taylor1, Lfunc_Taylor(L)));

  arb_clear(taylor1);
  Lfunc_clear(L);
}

// (3) supply after compute is fatal, by every route.
static void test_supply_after_compute(void) {
  double mus[] = {0, 1};

  // (3a) push after compute (void route: recorded, surfaced fatal at re-query
  // via the sticky supply_ecode; here we observe it on the next compute).
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
    assert(!fatal_error(ec));
    ec |= Lfunc_use_all_lpolys(L, cb, NULL);
    assert(!fatal_error(ec));
    ec |= Lfunc_compute(L);
    assert(!fatal_error(ec));
    acb_poly_t f; acb_poly_init(f); cb(f, 101, 2, 100, NULL);
    Lfunc_use_lpoly(L, 101, f); // void: records ERR_LIFECYCLE
    acb_poly_clear(f);
    Lerror_t e = Lfunc_compute(L);
    assert((e & ERR_LIFECYCLE) && fatal_error(e));
    Lfunc_clear(L);
  }

  // (3b) callback after compute (returns the code directly).
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
    assert(!fatal_error(ec));
    ec |= Lfunc_use_all_lpolys(L, cb, NULL);
    assert(!fatal_error(ec));
    ec |= Lfunc_compute(L);
    assert(!fatal_error(ec));
    Lerror_t e = Lfunc_use_all_lpolys(L, cb, NULL);
    assert((e & ERR_LIFECYCLE) && fatal_error(e));
    Lfunc_clear(L);
  }
}

int main(void) {
  test_compute_no_supply();
  test_double_compute();
  test_supply_after_compute();
  printf("compute_lifecycle_test: all tests passed\n");
  return 0;
}
