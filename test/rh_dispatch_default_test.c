// Default-verifier dispatch test. Sym^2 of the elliptic curve 11.a (degree 3,
// conductor 121, self_dual=YES). With the Turing default, turing_check_RH
// cannot confirm completeness at degree 3 and sets ERR_RH_ERROR (a warning).
// With the Buthe default, S = Wf + Winf - Ws lies in [0,0.98) and RH is
// confirmed, so ERR_RH_ERROR must be ABSENT. We assert on the flag and on a
// certified value (the first zero of L(Sym^2 11.a), LMFDB), never on stdout.
//
// FAILS on the pre-dispatch (Turing-default) tree, PASSES once Buthe is the
// default. Exit 0 = pass (asserts abort otherwise; the explicit `return 1`
// guards keep it failing under -DNDEBUG too).
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

static long g_ainvs[5];
static int g_sym;
static long g_bad, g_apbad;

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
  if (g_sym == 2) {
    long t = a * a - P;
    acb_poly_set_coeff_si(poly, 0, 1); acb_poly_set_coeff_si(poly, 1, -t);
    acb_poly_set_coeff_si(poly, 2, P * t); acb_poly_set_coeff_si(poly, 3, -(P * P * P));
  } else {
    acb_poly_set_coeff_si(poly, 0, 1); acb_poly_set_coeff_si(poly, 1, -a);
    acb_poly_set_coeff_si(poly, 2, P);
  }
}

int main(void) {
  long a11[5] = {0, -1, 1, -10, -20};
  double m3[3] = {0, 1, 0};
  for (int i = 0; i < 5; i++) g_ainvs[i] = a11[i];
  g_bad = 11; g_apbad = 1; g_sym = 2;

  Lerror_t ec = 0;
  char cd[] = ".";
  Lparams_t Lp = {0};
  Lp.degree = 3; Lp.conductor = 121; Lp.normalisation = 1.0; Lp.mus = m3;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.self_dual = YES; Lp.rank = DK; Lp.cache_dir = cd;

  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }

  // NO Lfunc_set_rh_method call: this exercises the DEFAULT verifier.
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("default-dispatch ecode = 0x%lx\n", (unsigned long)ec);
  if (ec & ERR_RH_ERROR) {
    fprintf(stderr, "FAIL: default verifier left ERR_RH_ERROR set for Sym^2(11.a)\n");
    Lfunc_clear(L);
    return 1;
  }
  assert((ec & ERR_RH_ERROR) == 0);

  // certify the first zero of L(Sym^2 11.a): t1 ~ 3.899... (empirically verified)
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  arb_t ref; arb_init(ref);
  arb_set_str(ref, "3.89928149477134478", 100);
  arb_add_error_2exp_si(ref, -30); // generous; the computed ball is far tighter
  int zok = arb_overlaps(zeros + 0, ref);
  arb_clear(ref);
  if (!zok) {
    fprintf(stderr, "FAIL: first zero did not overlap the LMFDB reference\n");
    Lfunc_clear(L);
    return 1;
  }
  assert(zok);

  Lfunc_clear(L);
  printf("PASS: default verifier is Buthe; RH confirmed and first zero certified.\n");
  return 0;
}
