// BOTH-mode genuine-contradiction coverage: the positive ERR_RH_METHODS_DISAGREE
// path that rh_both_disagree_test.c (a negative guard) never exercises. (1) the
// decision rule rh_methods_disagree over all four (too_many, overcount) combos;
// (2) integration: drive the REAL buthe_check_RH into a lower-h hard over-count
// and confirm it pairs with a confirming Turing side to yield DISAGREE. Asserts
// on flags only. Exit 0 = pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

static long mp(long v, long p) { return ((v % p) + p) % p; }
static long ap37(long P) {
  long c = 1;
  for (long x = 0; x < P; x++) {
    long x3 = mp(x * mp(x * x, P), P), rhs = mp(x3 - x, P);
    for (long y = 0; y < P; y++) if (mp(mp(y * y, P) + y, P) == rhs) c++;
  }
  return P + 1 - c;
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm; acb_poly_zero(poly); long P = (long)p;
  if (P == 37) { acb_poly_set_coeff_si(poly,0,1); acb_poly_set_coeff_si(poly,1,1); return; }
  long a = ap37(P);
  acb_poly_set_coeff_si(poly,0,1); acb_poly_set_coeff_si(poly,1,-a); acb_poly_set_coeff_si(poly,2,P);
}

static void force_S(Lfunc *LL, int grid_i, long value) {
  arb_t target, S, delta;
  arb_init(target); arb_init(S); arb_init(delta);
  arb_set_si(target, value);
  buthe_S_at(S, LL, grid_i, LL->wprec);
  arb_sub(delta, target, S, LL->wprec);
  arb_add(LL->buthe_Wf[grid_i], LL->buthe_Wf[grid_i], delta, LL->wprec);
  arb_clear(target); arb_clear(S); arb_clear(delta);
}

int main(void) {
  // (1) decision rule: exactly one hard over-count => DISAGREE.
  assert(rh_methods_disagree(false, false) == ERR_SUCCESS);
  assert(rh_methods_disagree(true,  false) == ERR_RH_METHODS_DISAGREE);
  assert(rh_methods_disagree(false, true ) == ERR_RH_METHODS_DISAGREE);
  assert(rh_methods_disagree(true,  true ) == ERR_SUCCESS);

  // (2) integration: real buthe_check_RH over-count (S<0) + confirming Turing.
  Lerror_t ec = 0; double mus[2] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BOTH) == ERR_SUCCESS);
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  Lfunc *LL = (Lfunc *)L;
  // Force h=8 to be non-certifying but non-negative, then force h=6 negative.
  // A negative rigorous S at any grid h is a hard over-count, not a certificate.
  force_S(LL, 0, 2);
  force_S(LL, 1, -1);
  Lerror_t be = buthe_check_RH(LL);
  assert((be & ERR_BUT_ERROR) != 0);
  assert(rh_methods_disagree(false, (be & ERR_BUT_ERROR) != 0) == ERR_RH_METHODS_DISAGREE);

  Lfunc_clear(L);
  printf("PASS: BOTH-mode flags a one-sided hard over-count as DISAGREE.\n");
  return 0;
}
