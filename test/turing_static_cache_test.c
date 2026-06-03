// Regression test for the turing.c static-cache bug.
//
// Invariant: an L-function's Turing RH verdict must not depend on what other
// L-functions were processed earlier in the same process. The degree-2
// elliptic curve 37.a1 (rank 1) confirms RH cleanly (ecode == 0) when it is the
// only object in a process. Here we process a DIFFERENT-degree object FIRST
// (Sym^2 of 11.a, degree 3) so that the static Turing caches in turing.c are
// initialised with degree-3 constants (h=B/16, t0=B/8, rlogpi, t0hbit with
// B=512/3), then compute 37.a1. On the broken library 37.a1 reuses those stale
// degree-3 constants (correct would be B=512/2) and spuriously returns
// ERR_RH_ERROR. The fix recomputes the constants per L-function.
//
// NOTE: the FIRST object to reach turing_check_RH must be the degree-3 one;
// if 37.a1 ran first it would itself seed the (then-correct) degree-2 cache and
// the bug would be masked. So this test deliberately orders sym2 before 37.a1.
//
// FAILS (assert) on the pre-fix library, PASSES on the post-fix library.

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

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
static Lerror_t run(long av[5], uint64_t cond, long bad, long apbad, int sym,
                    uint64_t deg, double norm, double *mus) {
  for (int i = 0; i < 5; i++) g_ainvs[i] = av[i];
  g_bad = bad; g_apbad = apbad; g_sym = sym;
  Lerror_t ec = 0; char cd[] = ".";
  Lparams_t Lp = {0};
  Lp.degree = deg; Lp.conductor = cond; Lp.normalisation = norm; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.self_dual = YES; Lp.rank = DK; Lp.cache_dir = cd;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  ec |= Lfunc_compute(L);
  Lfunc_clear(L);
  return ec;
}

int main(void) {
  long a11[5] = {0, -1, 1, -10, -20}, a37[5] = {0, 0, 1, -1, 0};
  double m2[2] = {0, 1}, m3[3] = {0, 1, 0};

  // Degree-3 object FIRST: seeds the static Turing caches with degree-3 values.
  (void)run(a11, 121, 11, 1, 2, 3, 1.0, m3); // Sym^2 of 11.a, degree 3

  // Then the degree-2 curve 37.a1. Run alone it confirms RH with ecode == 0.
  Lerror_t after = run(a37, 37, 37, -1, 1, 2, 0.5, m2);
  printf("37.a1 after deg-3 sym2: ecode=0x%lx (want 0x0)\n", (unsigned long)after);

  // The bug: `after` gains ERR_RH_ERROR on the broken library because the stale
  // degree-3 Turing constants are applied to this degree-2 zero count.
  assert((after & ERR_RH_ERROR) == 0);
  assert(after == 0);

  printf("PASS: 37.a1 RH verdict is order-independent.\n");
  return 0;
}
