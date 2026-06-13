/*
  Regression test: Lfunc_nmax signals failure with a 0 return.

  Lfunc_nmax returns 0 only when it cannot derive the coefficient bound M (a
  fatal internal error such as ERR_M_ERROR / ERR_OOM); a healthy object always
  returns nmax >= 1. An nmax-only caller (one that queries the bound without
  immediately supplying and computing, e.g. test/highdeg_check.cpp's nmax-query
  mode sizing a factor request) must treat a 0 return as an error, not as "no
  primes needed". This pins that contract: a pathological conductor must yield
  exactly 0 (a distinguishable error signal), while a healthy one yields >= 1.

  Through the public glfunc.h API only; exits 0 on success, asserts on failure.
*/
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

// A trivial well-formed degree-2 factor, so a supply call is not itself rejected
// for a bad argument; we only want to observe whether the nmax failure persists.
static void cb_one(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param) {
  (void)p; (void)d; (void)prec; (void)param;
  acb_poly_zero(poly);
  acb_poly_set_coeff_si(poly, 0, 1); // L_p(T) = 1: a valid (trivial) factor
}

int main(void) {
  double mus[] = {0, 1};

  // Healthy small object: nmax is a genuine positive bound.
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 11, 0.5, mus, &ec);
    assert(!fatal_error(ec) && L != NULL);
    uint64_t nmax = Lfunc_nmax(L);
    assert(nmax >= 1); // not a failure signal
    Lfunc_clear(L);
  }

  // Pathological conductor (2^60, degree 2): deriving M overflows the
  // coefficient array, so Lfunc_nmax cannot produce a bound and must return 0
  // (the failure signal), not a positive garbage value.
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, (uint64_t)1 << 60, 0.5, mus, &ec);
    assert(!fatal_error(ec) && L != NULL); // init itself succeeds
    uint64_t nmax = Lfunc_nmax(L);
    assert(nmax == 0); // failure is signalled as 0, distinguishable from a real bound

    // The failure is recorded on the object, so a caller that ignored the 0
    // return and proceeded to supply is still stopped fatally (here via a
    // callback route, which re-checks the recorded error and returns it).
    Lerror_t e = Lfunc_use_all_lpolys(L, cb_one, NULL);
    assert(fatal_error(e));
    Lfunc_clear(L);
  }

  printf("nmax_failure_test: all tests passed\n");
  return 0;
}
