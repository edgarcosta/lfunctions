// White-box test for the full Buthe RH-verification path (src/buthe.c), with
// the dynamic archimedean integral src/buthe_winf.c standing in for the old
// static gp/buthe_ints.out table.
//
// We build Sym^2 of the elliptic curve 11.a (degree 3, conductor 121,
// self_dual=YES) through the public API and run Lfunc_compute, so the active
// Turing verifier finds the zeros. We then reach into the internals
// (glfunc_internals.h), recover the Lfunc*, and call buthe_check_RH directly.
//
// For this degree-3 object Buthe's inequality is satisfied (the test sum
// S = Wf + Winf - Ws lies in [0, 0.98)), so buthe_check_RH must return
// ERR_SUCCESS: NO ERR_BUT_ERROR (S >= 0) and NO ERR_RH_ERROR (S < 0.98).
// (Turing, by contrast, miscounts for degree >= 3 and rejects this object;
// confirming Buthe accepts it is exactly the point of exercising this path.)
//
// The construction (point-counting good factors + a hardcoded bad factor,
// Sym^2 expansion of a_p) is copied from test/turing_static_cache_test.c.
//
// Assertions are on the Lerror_t flags returned by buthe_check_RH, never on
// stdout. Failure exits nonzero (the `if (... ) return 1;` guards) so the test
// still fails under -DNDEBUG where assert() is compiled out; the asserts are
// kept as belt-and-suspenders for assert-enabled builds.

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

static long g_ainvs[5];
static long g_bad, g_apbad;
static int g_sym;

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
  // Sym^2 of 11.a: a-invariants of 11.a1, conductor 121, bad prime 11
  // (a_p = 1 there for Sym^2), degree 3, normalisation 1.0, mus {0,1,0}.
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
  if (fatal_error(ec)) {
    fprintf(stderr, "Lfunc_init_advanced failed: ");
    fprint_errors(stderr, ec);
    return 2;
  }
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) {
    fprintf(stderr, "Lfunc_use_all_lpolys failed: ");
    fprint_errors(stderr, ec);
    return 2;
  }

  // Lfunc_compute runs the active (Turing) verifier and finds the zeros. For a
  // degree-3 object Turing miscounts and raises ERR_RH_ERROR as a *warning*;
  // that is not fatal, and we discard the verifier's verdict here -- we only
  // need the computed zeros so buthe_check_RH has something to sum over.
  ec |= Lfunc_compute(L);
  if (fatal_error(ec)) {
    fprintf(stderr, "Lfunc_compute failed (fatal): ");
    fprint_errors(stderr, ec);
    Lfunc_clear(L);
    return 2;
  }

  // White-box: run Buthe's verifier directly on the computed zeros.
  Lerror_t buthe = buthe_check_RH((Lfunc *) L);

  printf("buthe_check_RH returned 0x%lx\n", (unsigned long) buthe);
  if (buthe & ERR_BUT_ERROR)
    printf("  ERR_BUT_ERROR set: Buthe sum S < 0 (RH not confirmed).\n");
  if (buthe & ERR_RH_ERROR)
    printf("  ERR_RH_ERROR set: Buthe sum S >= 0.98 (missed zeros).\n");

  int ok = 1;
  if (buthe & ERR_BUT_ERROR) {
    fprintf(stderr, "FAIL: Buthe reported ERR_BUT_ERROR for Sym^2(11.a).\n");
    ok = 0;
  }
  if (buthe & ERR_RH_ERROR) {
    fprintf(stderr, "FAIL: Buthe reported ERR_RH_ERROR for Sym^2(11.a).\n");
    ok = 0;
  }
  if (buthe != ERR_SUCCESS) {
    fprintf(stderr, "FAIL: buthe_check_RH did not return ERR_SUCCESS (0x%lx).\n",
            (unsigned long) buthe);
    ok = 0;
  }

  assert((buthe & ERR_BUT_ERROR) == 0);
  assert((buthe & ERR_RH_ERROR) == 0);
  assert(buthe == ERR_SUCCESS);

  Lfunc_clear(L);

  if (!ok) {
    printf("buthe_path_test: FAILED\n");
    return 1;
  }
  printf("PASS: Buthe confirms RH (S<0.98) for Sym^2(11.a) via the dynamic W_infinity path.\n");
  return 0;
}
