// BOTH-mode false-disagreement guard. Sym^2 of 11.a (deg 3, self_dual=YES):
// Buthe CONFIRMS (S<0.98) but Turing cannot confirm at deg 3 and reports
// "too few" (NOT a hard over-count). Per the design, a verifier merely failing
// to certify is NOT a disagreement, so ERR_RH_METHODS_DISAGREE must be ABSENT
// and BOTH must return Buthe's (confirming) verdict: ERR_RH_ERROR absent too.
// FAILS on a naive disagree rule, PASSES on the genuine-contradiction rule.
// Asserts on flags only. Exit 0 = pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

static long g_ainvs[5];
static int g_sym; static long g_bad, g_apbad;
static long mod_p(long v, long p) { return ((v % p) + p) % p; }
static long ap_good(long p) {
  long a1 = g_ainvs[0], a2 = g_ainvs[1], a3 = g_ainvs[2], a4 = g_ainvs[3], a6 = g_ainvs[4];
  long c = 1;
  for (long x = 0; x < p; x++) {
    long x2 = mod_p(x * x, p), r = mod_p(x2 * x, p);
    r = mod_p(r + a2 * x2, p); r = mod_p(r + a4 * x, p); r = mod_p(r + a6, p);
    for (long y = 0; y < p; y++) if (mod_p(y * y + a1 * x * y + a3 * y, p) == r) c++;
  }
  return p + 1 - c;
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm; acb_poly_zero(poly); long P = (long)p;
  if (P == g_bad) {
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, (g_sym == 2) ? -1 : -g_apbad);
    return;
  }
  long a = ap_good(P);
  long t = a * a - P;
  acb_poly_set_coeff_si(poly, 0, 1); acb_poly_set_coeff_si(poly, 1, -t);
  acb_poly_set_coeff_si(poly, 2, P * t); acb_poly_set_coeff_si(poly, 3, -(P * P * P));
}

int main(void) {
  long a11[5] = {0, -1, 1, -10, -20};
  double m3[3] = {0, 1, 0};
  for (int i = 0; i < 5; i++) g_ainvs[i] = a11[i];
  g_bad = 11; g_apbad = 1; g_sym = 2;

  Lerror_t ec = 0; char cd[] = ".";
  Lparams_t Lp = {0};
  Lp.degree = 3; Lp.conductor = 121; Lp.normalisation = 1.0; Lp.mus = m3;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.self_dual = YES; Lp.rank = DK; Lp.cache_dir = cd;

  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  assert(Lfunc_set_rh_method(L, LFUNC_RH_BOTH) == ERR_SUCCESS);
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("BOTH(deg3) ecode = 0x%lx\n", (unsigned long)ec);
  int ok = 1;
  if (ec & ERR_RH_METHODS_DISAGREE) { fprintf(stderr,
      "FAIL: spurious ERR_RH_METHODS_DISAGREE (Turing too-few is not a disagreement)\n"); ok = 0; }
  if (ec & ERR_RH_ERROR) { fprintf(stderr,
      "FAIL: ERR_RH_ERROR set; BOTH should return Buthe's confirming verdict\n"); ok = 0; }
  assert((ec & ERR_RH_METHODS_DISAGREE) == 0);
  assert((ec & ERR_RH_ERROR) == 0);

  Lfunc_clear(L);
  if (!ok) return 1;
  printf("PASS: BOTH mode does not flag Turing's too-few as a disagreement.\n");
  return 0;
}
