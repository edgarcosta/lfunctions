// Buthe witness (public API): elliptic curve 37.a1 (deg 2, cond 37, self-dual,
// rank 1). With Buthe selected, RH is confirmed (no ERR_RH_ERROR/ERR_BUT_ERROR)
// and the numeric core is correct: |eps|=1 and the first zero overlaps the
// LMFDB value 5.0031700140066586953. Asserts on certified balls only. Exit 0 =
// pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

// 37.a1: ainvs [0,0,1,-1,0]. Bad prime 37 has local factor 1 + T (this is the
// LMFDB factor {37,{1,1}}; in the 1 - a_p T convention that is a_37 = -1).
// Good-prime local factor: 1 - a_p T + p T^2 (coeff[1] = -a_p).
// Use point-counting to derive a_p so signs are correct by construction.
static long mod_p37(long v, long p) { return ((v % p) + p) % p; }
static long ap_37(long P) {
  // ainvs [a1,a2,a3,a4,a6] = [0,0,1,-1,0]
  long c = 1;
  for (long x = 0; x < P; x++) {
    long x3 = mod_p37(x * mod_p37(x * x, P), P);
    long rhs = mod_p37(x3 - x + 0, P); // x^3 - x (a4=-1, a2=a6=0)
    for (long y = 0; y < P; y++)
      if (mod_p37(mod_p37(y * y, P) + y, P) == rhs) c++; // y^2 + a3*y = y^2+y
  }
  return P + 1 - c;
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm;
  acb_poly_zero(poly);
  long P = (long)p;
  if (P == 37) { // bad prime: local factor 1 + T (LMFDB {1,1}; a_37 = -1)
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, 1);
    return;
  }
  long a = ap_37(P);
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_poly_set_coeff_si(poly, 1, -a); // 1 - a_p T + p T^2
  acb_poly_set_coeff_si(poly, 2, P);
}

int main(void) {
  Lerror_t ec = 0;
  double mus[2] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BUTHE) == ERR_SUCCESS); // explicit

  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("witness ecode = 0x%lx\n", (unsigned long)ec);
  int ok = 1;
  if (ec & ERR_RH_ERROR) { fprintf(stderr, "FAIL: ERR_RH_ERROR (Buthe should confirm 37.a1)\n"); ok = 0; }
  if (ec & ERR_BUT_ERROR) { fprintf(stderr, "FAIL: ERR_BUT_ERROR\n"); ok = 0; }
  assert((ec & ERR_RH_ERROR) == 0);
  assert((ec & ERR_BUT_ERROR) == 0);

  // |eps| = 1
  acb_srcptr eps = Lfunc_sign(L);
  arb_t mag, one; arb_init(mag); arb_init(one);
  acb_abs(mag, eps, 100); arb_set_ui(one, 1);
  assert(arb_overlaps(mag, one));
  arb_clear(mag); arb_clear(one);

  // first zero overlaps LMFDB 5.0031700140066586953 (empirically verified: S=0.0053, rank 1)
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  arb_t zref; arb_init(zref);
  arb_set_str(zref, "5.0031700140066586953", 200);
  arb_add_error_2exp_si(zref, -40);
  if (!arb_overlaps(zeros + 0, zref)) { fprintf(stderr, "FAIL: first zero mismatch\n"); ok = 0; }
  assert(arb_overlaps(zeros + 0, zref));
  arb_clear(zref);

  Lfunc_clear(L);
  if (!ok) return 1;
  printf("PASS: Buthe confirms RH for 37.a1 with certified |eps|=1 and first zero.\n");
  return 0;
}
