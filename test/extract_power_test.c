// Hermetic extraction test: builds E^2 and E^3 from an elliptic curve's
// Euler factors and checks extraction recovers E and assembles L = E^k.
// Asserts on certified balls only. Build with assertions on.
#include <assert.h>
#include <stdio.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

// 37.a1 Euler polynomials L_p(T) = 1 + c1*T + p*T^2, stored as {p, c1, p};
// the bad prime 37 is {1, 1} (handled in Ep below). These c1 are the POLYNOMIAL
// T^1 coefficients (c1 = -a_p, NOT the curve's a_p), copied VERBATIM from the
// euler_factors map in examples/ec_37.a1.cpp. Do not negate them.
static const long ap37[][3] = {
  {2,2,2},{3,3,3},{5,2,5},{7,1,7},{11,5,11},{13,2,13},{17,0,17},
  {19,0,19},{23,-2,23},{29,-6,29},{31,4,31},{41,9,41},{43,-2,43},{47,9,47},
  {53,-1,53},{59,-8,59},{61,8,61},{67,-8,67},{71,-9,71},{73,1,73},{79,-4,79},
  {83,15,83},{89,-4,89},{97,-4,97},{101,-3,101},{103,-18,103},{107,12,107},
  {109,16,109},{113,18,113},{127,-1,127},{131,12,131},{137,6,137},{139,-4,139},
  {149,5,149},  // one good prime > nmax(37.a)=142, for the out-of-order supply test
};
// Sanity check against examples/ec_37.a1.cpp: p=2 -> {1,2,2}, p=23 -> {1,-2,23}.
// p=149 a_p=-5 (T-coeff 5) computed by point count on y^2+y=x^3-x mod 149.

static int k_param; // 2 or 3: the callback raises E_p to this power

static void Ep(acb_poly_t poly, uint64_t p) // E_p(T) for 37.a1, degree 2 (or 1 at p=37)
{
  acb_poly_zero(poly);
  if (p == 37) { acb_poly_set_coeff_si(poly,0,1); acb_poly_set_coeff_si(poly,1,1); return; }
  for (size_t i = 0; i < sizeof(ap37)/sizeof(ap37[0]); i++)
    if ((uint64_t)ap37[i][0] == p) {
      acb_poly_set_coeff_si(poly,0,1);
      acb_poly_set_coeff_si(poly,1,ap37[i][1]);
      acb_poly_set_coeff_si(poly,2,ap37[i][2]);
      return;
    }
  // unknown prime: leave zero -> Lfunc_use_all_lpolys short-circuits (fine; see plan)
}

// callback for L = E^k: supply (E_p)^k
static void lk_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  (void)d; (void)prec; (void)param;
  acb_poly_t ep; acb_poly_init(ep);
  Ep(ep, p);
  if (acb_poly_is_zero(ep)) { acb_poly_zero(poly); acb_poly_clear(ep); return; }
  acb_poly_pow_ui(poly, ep, (ulong)k_param, prec);
  acb_poly_clear(ep);
}

// callback for E itself: supply E_p
static void e_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  (void)d; (void)prec; (void)param;
  Ep(poly, p);
}

// callback for I3: like lk_callback but makes one coefficient INEXACT (a tiny ball).
// The k-th-power certificate is rigorous only for exact factors, so an inexact one
// must be refused even though its root^k still overlaps the wide ball.
static void lk_inexact_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  lk_callback(poly, p, d, prec, param);
  if (!acb_poly_is_zero(poly))
    arb_add_error_2exp_si(acb_realref(acb_poly_get_coeff_ptr(poly, 0)), -150);
}

int main(void)
{
  Lerror_t ec = ERR_SUCCESS;
  double mus[2] = {0,1};

  // (Task 2) A primitive L (E itself) reports zero factors.
  Lfunc_t E = Lfunc_init(2, 37, 0.5, mus, &ec);
  assert(!fatal_error(ec));
  ec |= Lfunc_use_all_lpolys(E, e_callback, NULL);
  ec |= Lfunc_compute(E);
  assert(!fatal_error(ec));
  Lfunc_t *facs = NULL; uint64_t *mults = NULL;
  assert(Lfunc_factors(E, &facs, &mults) == 0);
  Lfunc_clear(E);

  printf("extract_power_test: Task 2 OK\n");

  // Build E^2 with extract_powers = NO -> ERR_POWER.
  k_param = 2;
  Lparams_t Lp2 = { .degree = 4, .conductor = 37u*37u, .normalisation = 0.5,
                    .mus = (double[]){0,0,1,1}, .target_prec = 100, .wprec = 0,
                    .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                    .extract_powers = NO };
  Lerror_t e_off = ERR_SUCCESS;
  Lfunc_t Loff = Lfunc_init_advanced(&Lp2, &e_off);
  assert(!fatal_error(e_off));
  e_off |= Lfunc_use_all_lpolys(Loff, lk_callback, NULL); // INSUFF_EULER warning is fine
  e_off |= Lfunc_compute(Loff);
  assert(e_off & ERR_POWER);          // opted out: rejected
  Lfunc_clear(Loff);

  // Reference: compute E directly.
  Lfunc_t Eref = Lfunc_init(2, 37, 0.5, mus, &ec);
  ec = ERR_SUCCESS;
  ec |= Lfunc_use_all_lpolys(Eref, e_callback, NULL);
  ec |= Lfunc_compute(Eref);
  assert(!fatal_error(ec));
  int64_t rankE = Lfunc_rank(Eref);
  arb_srcptr zerosE = Lfunc_zeros(Eref, 0);

  // Build E^2 with extract_powers = YES -> extract & assemble.
  Lp2.extract_powers = YES;
  Lerror_t e2 = ERR_SUCCESS;
  Lfunc_t L2 = Lfunc_init_advanced(&Lp2, &e2);
  assert(!fatal_error(e2));
  e2 |= Lfunc_use_all_lpolys(L2, lk_callback, NULL);
  e2 |= Lfunc_compute(L2);
  assert(!fatal_error(e2));

  // factor accessor returns (M, 2); M overlaps E
  Lfunc_t *f2 = NULL; uint64_t *m2 = NULL;
  assert(Lfunc_factors(L2, &f2, &m2) == 1);
  assert(m2[0] == 2);
  assert(Lfunc_rank(f2[0]) == rankE);
  // assembled rank, zeros (doubled), sign
  assert(Lfunc_rank(L2) == 2*rankE);
  arb_srcptr z2 = Lfunc_zeros(L2, 0);
  assert(arb_overlaps((arb_ptr)(z2+0), (arb_ptr)(zerosE+0)));
  assert(arb_overlaps((arb_ptr)(z2+1), (arb_ptr)(zerosE+0))); // repeated k=2 times
  assert(arb_overlaps((arb_ptr)(z2+2), (arb_ptr)(zerosE+1)));
  acb_t s2; acb_init(s2);
  acb_pow_ui(s2, Lfunc_sign(Eref), 2, 100);           // sign(E)^2
  assert(acb_overlaps(s2, Lfunc_sign(L2)));
  { arb_t am; arb_init(am); acb_abs(am, Lfunc_sign(L2), 100);
    arb_sub_ui(am, am, 1, 100); assert(arb_contains_zero(am)); arb_clear(am); } // |eps|=1
  acb_clear(s2);

  // (Task 4) Special value: L2(1.5) == Eref(1.5)^2
  acb_t vL, vE, vEk; acb_init(vL); acb_init(vE); acb_init(vEk);
  Lerror_t sv = ERR_SUCCESS;
  sv |= Lfunc_special_value(vE, Eref, 1.5, 0.0);
  sv |= Lfunc_special_value(vL, L2, 1.5, 0.0);
  assert(!fatal_error(sv));
  acb_pow_ui(vEk, vE, 2, 100);          // L(s) = E(s)^2
  assert(acb_overlaps(vL, vEk));
  acb_clear(vL); acb_clear(vE); acb_clear(vEk);

  // (Task / spec 13) assembled leading Taylor coefficient: L_d(L2) == L_d(Eref)^2
  { arb_t tk; arb_init(tk);
    arb_pow_ui(tk, Lfunc_Taylor(Eref), 2, 100);
    assert(arb_overlaps((arb_ptr)Lfunc_Taylor(L2), tk));
    arb_clear(tk); }

  Lfunc_clear(L2);
  Lfunc_clear(Eref);

  printf("extract_power_test: Task 3 OK\n");

  // ---- E^3 (degree 6, odd k) ----
  k_param = 3;
  Lfunc_t Eref3 = Lfunc_init(2, 37, 0.5, mus, &ec);
  ec = ERR_SUCCESS;
  ec |= Lfunc_use_all_lpolys(Eref3, e_callback, NULL);
  ec |= Lfunc_compute(Eref3);
  assert(!fatal_error(ec));

  Lparams_t Lp3 = { .degree = 6, .conductor = 50653, .normalisation = 0.5,
                    .mus = (double[]){0,0,0,1,1,1}, .target_prec = 100, .wprec = 0,
                    .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                    .extract_powers = YES };
  Lerror_t e3 = ERR_SUCCESS;
  Lfunc_t L3 = Lfunc_init_advanced(&Lp3, &e3);
  assert(!fatal_error(e3));
  e3 |= Lfunc_use_all_lpolys(L3, lk_callback, NULL);
  e3 |= Lfunc_compute(L3);
  assert(!fatal_error(e3));

  Lfunc_t *f3 = NULL; uint64_t *m3 = NULL;
  assert(Lfunc_factors(L3, &f3, &m3) == 1 && m3[0] == 3);
  assert(Lfunc_rank(L3) == 3*Lfunc_rank(Eref3));
  arb_srcptr zE3 = Lfunc_zeros(Eref3, 0);
  arb_srcptr z3  = Lfunc_zeros(L3, 0);
  assert(arb_overlaps((arb_ptr)(z3+0), (arb_ptr)(zE3+0)));
  assert(arb_overlaps((arb_ptr)(z3+1), (arb_ptr)(zE3+0)));
  assert(arb_overlaps((arb_ptr)(z3+2), (arb_ptr)(zE3+0))); // tripled
  assert(arb_overlaps((arb_ptr)(z3+3), (arb_ptr)(zE3+1)));
  acb_t s3; acb_init(s3); acb_pow_ui(s3, Lfunc_sign(Eref3), 3, 100);
  assert(acb_overlaps(s3, Lfunc_sign(L3))); acb_clear(s3); // sign(E)^3 = -1
  acb_t vL3, vE3, vE3k; acb_init(vL3); acb_init(vE3); acb_init(vE3k);
  Lerror_t sv3 = Lfunc_special_value(vE3, Eref3, 1.5, 0.0)
               | Lfunc_special_value(vL3, L3, 1.5, 0.0);
  assert(!fatal_error(sv3));
  acb_pow_ui(vE3k, vE3, 3, 100);
  assert(acb_overlaps(vL3, vE3k));
  acb_clear(vL3); acb_clear(vE3); acb_clear(vE3k);

  // (Task / spec 13) assembled leading Taylor coefficient: L_d(L3) == L_d(Eref3)^3
  { arb_t tk; arb_init(tk);
    arb_pow_ui(tk, Lfunc_Taylor(Eref3), 3, 100);
    assert(arb_overlaps((arb_ptr)Lfunc_Taylor(L3), tk));
    arb_clear(tk); }

  Lfunc_clear(L3); Lfunc_clear(Eref3);

  printf("extract_power_test: Task 5 OK\n");

  // ---- out-of-order supply (regression for extract_and_assemble's prime loop) ----
  // The highdeg harness supplies good factors ascending, then appends the BAD
  // factors at the tail (e.g. p=2,7 after p up to ~1900 for 196.a). So L's retained
  // factors are NOT in prime order. extract_and_assemble must feed M every retained
  // prime <= nmax(M), regardless of position; a `break` at the first p > nmax(M)
  // would drop the trailing low primes and compute M (here 37.a, nmax 142) WITHOUT
  // its bad factor at 37 -> wrong zeros/Taylor. We reproduce that layout here:
  // supply E^2's good primes ascending (incl. p=149 > 142) via Lfunc_use_lpoly, then
  // the bad prime 37 LAST.
  k_param = 2;
  Lfunc_t ErefO = Lfunc_init(2, 37, 0.5, mus, &ec);
  ec = ERR_SUCCESS;
  ec |= Lfunc_use_all_lpolys(ErefO, e_callback, NULL);
  ec |= Lfunc_compute(ErefO);
  assert(!fatal_error(ec));
  arb_srcptr zerosEO = Lfunc_zeros(ErefO, 0);

  Lparams_t LpO = { .degree = 4, .conductor = 37u*37u, .normalisation = 0.5,
                    .mus = (double[]){0,0,1,1}, .target_prec = 100, .wprec = 0,
                    .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                    .extract_powers = YES };
  Lerror_t eO = ERR_SUCCESS;
  Lfunc_t LO = Lfunc_init_advanced(&LpO, &eO);
  assert(!fatal_error(eO));
  {
    acb_poly_t fp; acb_poly_init(fp);
    uint64_t nmaxO = Lfunc_nmax(LO);
    // good primes (everything in ap37 except the bad prime 37), ascending,
    // each raised to k=2; this includes p=149 > nmax(37.a)=142.
    for (size_t i = 0; i < sizeof(ap37)/sizeof(ap37[0]); i++) {
      uint64_t p = (uint64_t)ap37[i][0];
      if (p == 37 || p > nmaxO) continue;
      acb_poly_t ep; acb_poly_init(ep); Ep(ep, p);
      acb_poly_pow_ui(fp, ep, 2, 100);
      Lfunc_use_lpoly(LO, p, fp);
      acb_poly_clear(ep);
    }
    // bad prime 37 supplied LAST (out of order), (1+T)^2:
    if (37 <= nmaxO) {
      acb_poly_t ep; acb_poly_init(ep); Ep(ep, 37);
      acb_poly_pow_ui(fp, ep, 2, 100);
      Lfunc_use_lpoly(LO, 37, fp);
      acb_poly_clear(ep);
    }
    acb_poly_clear(fp);
  }
  eO |= Lfunc_compute(LO);
  assert(!fatal_error(eO));
  Lfunc_t *fO = NULL; uint64_t *mO = NULL;
  assert(Lfunc_factors(LO, &fO, &mO) == 1 && mO[0] == 2);
  // The assembled first zero must equal 37.a's first zero. With the `break` bug,
  // M is computed without its bad factor at 37 and this overlap fails.
  arb_srcptr zO = Lfunc_zeros(LO, 0);
  assert(arb_overlaps((arb_ptr)(zO+0), (arb_ptr)(zerosEO+0)));
  assert(arb_overlaps((arb_ptr)(zO+1), (arb_ptr)(zerosEO+0))); // doubled
  assert(arb_overlaps((arb_ptr)(zO+2), (arb_ptr)(zerosEO+1)));
  Lfunc_clear(LO); Lfunc_clear(ErefO);

  printf("extract_power_test: out-of-order supply OK\n");

  // ---- I2: non-k-divisible mus must be rejected (soundness gate) ----
  // A degree-4 L with mus={0,0,0,1} has sorted blocks [0,0] and [0,1] (k=2),
  // which are unequal: L's gamma factors disagree with any pure square, so even
  // if the Euler factors are perfect squares we must NOT extract. Supply the same
  // perfect-square factors as E^2 via lk_callback; expect ERR_POWER.
  k_param = 2;
  Lparams_t LpMu = { .degree = 4, .conductor = 37u*37u, .normalisation = 0.5,
                     .mus = (double[]){0,0,0,1}, .target_prec = 100, .wprec = 0,
                     .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                     .extract_powers = YES };
  Lerror_t eMu = ERR_SUCCESS;
  Lfunc_t LMu = Lfunc_init_advanced(&LpMu, &eMu);
  assert(!fatal_error(eMu));
  eMu |= Lfunc_use_all_lpolys(LMu, lk_callback, NULL); // INSUFF_EULER warning is fine
  eMu |= Lfunc_compute(LMu);
  assert(eMu & ERR_POWER);            // mus not k-divisible -> refuse to extract
  Lfunc_clear(LMu);

  printf("extract_power_test: I2 non-k-divisible mus OK\n");

  // ---- I3: inexact (wide-ball) Euler factors must not be certified (rigor) ----
  // E^2 with one coefficient carrying a tiny error ball. The certificate proves
  // f = M_p^k by ball containment, which is meaningful only for exact factors; a
  // wide-ball non-power can false-accept. Expect ERR_POWER (refuse to certify).
  k_param = 2;
  Lparams_t LpIn = { .degree = 4, .conductor = 37u*37u, .normalisation = 0.5,
                     .mus = (double[]){0,0,1,1}, .target_prec = 100, .wprec = 0,
                     .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                     .extract_powers = YES };
  Lerror_t eIn = ERR_SUCCESS;
  Lfunc_t LIn = Lfunc_init_advanced(&LpIn, &eIn);
  assert(!fatal_error(eIn));
  eIn |= Lfunc_use_all_lpolys(LIn, lk_inexact_callback, NULL); // INSUFF_EULER warning is fine
  eIn |= Lfunc_compute(LIn);
  assert(eIn & ERR_POWER);            // inexact factor -> refuse to certify
  Lfunc_clear(LIn);

  printf("extract_power_test: I3 inexact factor OK\n");

  return 0;
}
