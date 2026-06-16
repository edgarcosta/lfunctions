/*
  Regression test: Lfunc_reduce_nmax must reject a zero reduction.

  Lfunc_reduce_nmax(L, 0) used to be accepted (0 < M is true), setting the
  internal coefficient count M to 0. A later M_error divides by L->M (see
  src/error.c), so a subsequent Lfunc_compute divided by zero and produced a
  NaN/infinite error bound rather than a meaningful result. M = 1 (a single
  coefficient a_1 = 1) is the smallest defensible bound, so 0 must be rejected
  cleanly.

  We assert on the public-API return value (false) and the recorded fatal
  supply error, and on a positive boundary value (the largest accepted nmax,
  M - 1) still succeeding. Through the public glfunc.h API only; exits 0 on
  success, asserts on failure.
*/
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include "glfunc.h"

int main(void) {
  double mus[] = {0, 1};

  // Zero reduction is rejected: returns false and records a fatal supply error
  // that Lfunc_compute then refuses to run on.
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
    assert(!fatal_error(ec));
    uint64_t M = Lfunc_nmax(L);
    assert(M > 1);

    bool ok = Lfunc_reduce_nmax(L, 0);
    assert(!ok); // must NOT accept a zero reduction

    // The rejection is sticky: it is recorded on the object, so a caller that
    // ignores the false return still cannot compute a bogus result.
    Lerror_t e = Lfunc_compute(L);
    assert(fatal_error(e));
    Lfunc_clear(L);
  }

  // The positive boundary value M - 1 (the largest nmax strictly below M) is
  // still accepted, so the fix rejects only 0, not the whole low end.
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
    assert(!fatal_error(ec));
    uint64_t M = Lfunc_nmax(L);
    assert(M > 1);

    bool ok = Lfunc_reduce_nmax(L, M - 1);
    assert(ok); // a genuine, positive reduction is accepted
    Lfunc_clear(L);
  }

  // nmax >= M is still rejected (unchanged contract: cannot increase the bound).
  {
    Lerror_t ec = ERR_SUCCESS;
    Lfunc_t L = Lfunc_init(2, 35, 0.0, mus, &ec);
    assert(!fatal_error(ec));
    uint64_t M = Lfunc_nmax(L);
    assert(!Lfunc_reduce_nmax(L, M));
    Lfunc_clear(L);
  }

  printf("reduce_nmax_zero_test: all tests passed\n");
  return 0;
}
