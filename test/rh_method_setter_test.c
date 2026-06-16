// Unit test for the RH-method setter contract (Lfunc_set_rh_method) and the
// default verifier (LFUNC_RH_BUTHE). We only init an L-function (no compute),
// then exercise the setter's accept/reject/clamp behavior, and finally that the
// post-compute call is rejected. Assertions are on Lerror_t flags / the
// returned value, never on stdout. Exit 0 = pass (asserts abort otherwise).
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

// trivial degree-2 self-dual object so init is cheap; we never call compute
// except once at the very end to test the post-compute rejection.
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm; (void)p;
  acb_poly_zero(poly);
  acb_poly_set_coeff_si(poly, 0, 1); // 1 + 0*T + p*T^2 stand-in is fine for a
  acb_poly_set_coeff_si(poly, 2, (long)p); // self-dual deg-2 shape; values irrelevant here
}

int main(void) {
  Lerror_t ec = 0;
  double mus[2] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 11, 0.5, mus, &ec);
  assert(!fatal_error(ec));

  // default must be Buthe
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_BUTHE);

  // valid sets succeed and stick
  assert(Lfunc_set_rh_method(L, LFUNC_RH_TURING) == ERR_SUCCESS);
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_TURING);
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BOTH) == ERR_SUCCESS);
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_BOTH);
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BUTHE) == ERR_SUCCESS);
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_BUTHE);

  // an out-of-range method is rejected and leaves rh_method unchanged
  assert(Lfunc_set_rh_method(L, (Lfunc_rh_method)99) == ERR_BAD_RH_METHOD);
  assert(fatal_error(ERR_BAD_RH_METHOD));
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_BUTHE);

  // after compute, the setter is rejected
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  assert(!fatal_error(ec));
  ec |= Lfunc_compute(L); // verdict/warnings irrelevant; we only need computed=true
  assert(((Lfunc *)L)->computed);
  assert(Lfunc_set_rh_method(L, LFUNC_RH_TURING) == ERR_BAD_RH_METHOD);
  assert(((Lfunc *)L)->rh_method == LFUNC_RH_BUTHE); // unchanged

  Lfunc_clear(L);
  printf("PASS: rh_method default + setter contract.\n");
  return 0;
}
