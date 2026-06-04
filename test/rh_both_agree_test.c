// BOTH-mode agreement (deg 2). 37.a1: both Turing and Buthe confirm RH, so
// neither ERR_RH_ERROR nor ERR_RH_METHODS_DISAGREE is set, and the returned
// verdict (Buthe's) confirms. Asserts on flags + |eps|=1. Exit 0 = pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

// 37.a1: ainvs [0,0,1,-1,0]. Bad prime 37 has local factor 1 + T (LMFDB
// {37,{1,1}}; in the 1 - a_p T convention that is a_37 = -1).
// Good-prime local factor: 1 - a_p T + p T^2 (coeff[1] = -a_p).
// Use point-counting to derive a_p so signs are correct by construction.
static long mod_p37b(long v, long p) { return ((v % p) + p) % p; }
static long ap_37b(long P) {
  long c = 1;
  for (long x = 0; x < P; x++) {
    long x3 = mod_p37b(x * mod_p37b(x * x, P), P);
    long rhs = mod_p37b(x3 - x + 0, P);
    for (long y = 0; y < P; y++)
      if (mod_p37b(mod_p37b(y * y, P) + y, P) == rhs) c++;
  }
  return P + 1 - c;
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm;
  acb_poly_zero(poly);
  long P = (long)p;
  if (P == 37) {
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, 1); // 1 + T (LMFDB {1,1}; a_37 = -1)
    return;
  }
  long a = ap_37b(P);
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_poly_set_coeff_si(poly, 1, -a); // 1 - a_p T + p T^2
  acb_poly_set_coeff_si(poly, 2, P);
}

int main(void) {
  Lerror_t ec = 0;
  double mus[2] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BOTH) == ERR_SUCCESS);

  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("BOTH(deg2) ecode = 0x%lx\n", (unsigned long)ec);
  int ok = 1;
  if (ec & ERR_RH_METHODS_DISAGREE) { fprintf(stderr, "FAIL: spurious disagree flag\n"); ok = 0; }
  if (ec & ERR_RH_ERROR) { fprintf(stderr, "FAIL: ERR_RH_ERROR (both should confirm)\n"); ok = 0; }
  assert((ec & ERR_RH_METHODS_DISAGREE) == 0);
  assert((ec & ERR_RH_ERROR) == 0);

  acb_srcptr eps = Lfunc_sign(L);
  arb_t mag, one; arb_init(mag); arb_init(one);
  acb_abs(mag, eps, 100); arb_set_ui(one, 1);
  assert(arb_overlaps(mag, one));
  arb_clear(mag); arb_clear(one);

  Lfunc_clear(L);
  if (!ok) return 1;
  printf("PASS: BOTH mode -- Turing and Buthe agree, RH confirmed for 37.a1.\n");
  return 0;
}
