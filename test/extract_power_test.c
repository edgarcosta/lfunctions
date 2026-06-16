// Hermetic extraction test: builds E^2 and E^3 from an elliptic curve's
// Euler factors and checks extraction recovers E and assembles L = E^k.
// Asserts on certified balls only. Build with assertions on.
#include <assert.h>
#include <stdio.h>
#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/fmpz_poly.h>
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

static void Ep_z(fmpz_poly_t poly, uint64_t p) // exact E_p(T)
{
  fmpz_poly_zero(poly);
  if (p == 37) { fmpz_poly_set_coeff_ui(poly,0,1); fmpz_poly_set_coeff_ui(poly,1,1); return; }
  for (size_t i = 0; i < sizeof(ap37)/sizeof(ap37[0]); i++)
    if ((uint64_t)ap37[i][0] == p) {
      fmpz_poly_set_coeff_ui(poly,0,1);
      fmpz_poly_set_coeff_si(poly,1,ap37[i][1]);
      fmpz_poly_set_coeff_ui(poly,2,(ulong)ap37[i][2]);
      return;
    }
  // unknown prime: leave zero -> Lfunc_use_all_lpolys_fmpz short-circuits
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

// exact callback for L = E^k: supply (E_p)^k through fmpz provenance
static void lk_callback_z(fmpz_poly_t poly, uint64_t p, int d, void *param)
{
  (void)d; (void)param;
  fmpz_poly_t ep; fmpz_poly_init(ep);
  Ep_z(ep, p);
  if (fmpz_poly_is_zero(ep)) { fmpz_poly_zero(poly); fmpz_poly_clear(ep); return; }
  fmpz_poly_pow(poly, ep, (ulong)k_param);
  fmpz_poly_clear(ep);
}

static void lk_first_prime_only_callback_z(fmpz_poly_t poly, uint64_t p, int d, void *param)
{
  if (p == 2)
    lk_callback_z(poly, p, d, param);
  else
    fmpz_poly_zero(poly);
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

static void lk_badprime_callback_z(fmpz_poly_t poly, uint64_t p, int d, void *param)
{
  lk_callback_z(poly, p, d, param);
  if (p == 37) {
    fmpz_poly_zero(poly);
    fmpz_poly_set_coeff_ui(poly, 0, 1);
    fmpz_poly_set_coeff_si(poly, 1, -3);
    fmpz_poly_set_coeff_ui(poly, 2, 1);
  }
}

// Conductor-1 exact square with zero T^1 coefficient: L_p = (1 + T^2)^2.
// The power guard's moment signal is zero here, so extraction depends on treating
// conductor 1 as an exact k-th-power conductor and then certifying the fmpz roots.
static void conductor1_square_callback_z(fmpz_poly_t poly, uint64_t p, int d, void *param)
{
  (void)p; (void)d; (void)param;
  fmpz_poly_t root;
  fmpz_poly_init(root);
  fmpz_poly_set_coeff_ui(root, 0, 1);
  fmpz_poly_set_coeff_ui(root, 2, 1);
  fmpz_poly_pow(poly, root, 2);
  fmpz_poly_clear(root);
}

// callback for the complex (non-self-dual) boundary: build M with a COMPLEX linear
// coefficient (i times E_p's), so L = M^2 has complex Euler factors but the same 2nd
// moment and perfect-square conductor as E^2 (detection still finds k=2). The exact
// integer certificate requires real-integer coefficients, so it must refuse this with
// ERR_POWER -- extraction supports rational-integer L-functions only (complex
// extraction is a tracked follow-up).
static void lk_complex_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  (void)d; (void)param;
  acb_poly_t mp; acb_poly_init(mp);
  Ep(mp, p);
  if (acb_poly_is_zero(mp)) { acb_poly_zero(poly); acb_poly_clear(mp); return; }
  acb_t c1; acb_init(c1);
  acb_poly_get_coeff_acb(c1, mp, 1);
  acb_mul_onei(c1, c1);                 // i * (linear coeff): M is non-self-dual (complex a_p)
  acb_poly_set_coeff_acb(mp, 1, c1);
  acb_clear(c1);
  acb_poly_pow_ui(poly, mp, 2, prec);   // L_p = M_p^2 (complex)
  acb_poly_clear(mp);
}

static Lerror_t use_power_factor_fmpz(Lfunc_t L, uint64_t p, ulong k)
{
  fmpz_poly_t ep, fp;
  fmpz_poly_init(ep);
  fmpz_poly_init(fp);
  Ep_z(ep, p);
  fmpz_poly_pow(fp, ep, k);
  Lerror_t e = Lfunc_use_lpoly_fmpz(L, p, fp);
  fmpz_poly_clear(fp);
  fmpz_poly_clear(ep);
  return e;
}

static void use_power_factor_acb(Lfunc_t L, uint64_t p, ulong k)
{
  acb_poly_t ep, fp;
  acb_poly_init(ep);
  acb_poly_init(fp);
  Ep(ep, p);
  acb_poly_pow_ui(fp, ep, k, 100);
  Lfunc_use_lpoly(L, p, fp);
  acb_poly_clear(fp);
  acb_poly_clear(ep);
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
  e2 |= Lfunc_use_all_lpolys_fmpz(L2, lk_callback_z, NULL);
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

  // Exact-looking factors supplied only through the acb API are not exact
  // provenance for extraction. The guard still detects the square, but
  // power_extract_prepare must refuse to certify it.
  Lp2.extract_powers = YES;
  Lerror_t eA = ERR_SUCCESS;
  Lfunc_t LA = Lfunc_init_advanced(&Lp2, &eA);
  assert(!fatal_error(eA));
  eA |= Lfunc_use_all_lpolys(LA, lk_callback, NULL);
  eA |= Lfunc_compute(LA);
  assert(eA & ERR_POWER);
  Lfunc_t *fA = NULL; uint64_t *mA = NULL;
  assert(Lfunc_factors(LA, &fA, &mA) == 0);
	  Lfunc_clear(LA);

	  printf("extract_power_test: acb-only extraction refused OK\n");

	  // fmpz callback zero sentinel: first missing exact factor stops supply and
	  // reports ERR_INSUFF_EULER.
	  {
	    Lerror_t eS = ERR_SUCCESS;
	    Lfunc_t LS = Lfunc_init_advanced(&Lp2, &eS);
	    assert(!fatal_error(eS));
	    eS |= Lfunc_use_all_lpolys_fmpz(LS, lk_first_prime_only_callback_z, NULL);
	    assert(eS & ERR_INSUFF_EULER);
	    assert(!fatal_error(eS));
	    eS |= Lfunc_compute(LS);
	    assert(eS & ERR_INSUFF_EULER);
	    Lfunc_clear(LS);
	  }

	  printf("extract_power_test: fmpz callback sentinel OK\n");

	  // Mixed provenance must not extract: every retained factor must have exact
	  // fmpz provenance, regardless of supply order.
	  {
	    Lerror_t eMix1 = ERR_SUCCESS;
	    Lfunc_t LMix1 = Lfunc_init_advanced(&Lp2, &eMix1);
	    assert(!fatal_error(eMix1));
	    (void)Lfunc_nmax(LMix1);
	    eMix1 |= use_power_factor_fmpz(LMix1, 2, 2);
	    use_power_factor_acb(LMix1, 3, 2);
	    eMix1 |= Lfunc_compute(LMix1);
	    assert(eMix1 & ERR_POWER);
	    Lfunc_t *fm = NULL; uint64_t *mm = NULL;
	    assert(Lfunc_factors(LMix1, &fm, &mm) == 0);
	    Lfunc_clear(LMix1);

	    Lerror_t eMix2 = ERR_SUCCESS;
	    Lfunc_t LMix2 = Lfunc_init_advanced(&Lp2, &eMix2);
	    assert(!fatal_error(eMix2));
	    (void)Lfunc_nmax(LMix2);
	    use_power_factor_acb(LMix2, 2, 2);
	    eMix2 |= use_power_factor_fmpz(LMix2, 3, 2);
	    eMix2 |= Lfunc_compute(LMix2);
	    assert(eMix2 & ERR_POWER);
	    assert(Lfunc_factors(LMix2, &fm, &mm) == 0);
	    Lfunc_clear(LMix2);
	  }

  printf("extract_power_test: mixed provenance refused OK\n");

  // Conductor 1 is still an exact square conductor: exact supplied powers must
  // enter extraction instead of falling through as an ordinary primitive object.
  // A previous temporary reproducer for this path returned ERR_POWER=0 and
  // n_factors=0; the factor assertion below catches that bypass directly.
  {
    Lparams_t LpOne = { .degree = 4, .conductor = 1, .normalisation = 0.0,
                        .mus = (double[]){0,0,1,1}, .target_prec = 100,
                        .wprec = 0, .gprec = 0, .self_dual = YES, .rank = DK,
                        .cache_dir = ".", .extract_powers = YES };
    Lerror_t eOne = ERR_SUCCESS;
    Lfunc_t LOne = Lfunc_init_advanced(&LpOne, &eOne);
    assert(!fatal_error(eOne));
    eOne |= Lfunc_use_all_lpolys_fmpz(LOne, conductor1_square_callback_z, NULL);
    eOne |= Lfunc_compute(LOne);
    assert(!fatal_error(eOne));
    Lfunc_t *fOne = NULL; uint64_t *mOne = NULL;
    assert(Lfunc_factors(LOne, &fOne, &mOne) == 1);
    assert(mOne[0] == 2);
    assert(Lfunc_rank(LOne) == 2*Lfunc_rank(fOne[0]));
    Lfunc_clear(LOne);
  }

  printf("extract_power_test: conductor-1 exact square extracted OK\n");

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
  e3 |= Lfunc_use_all_lpolys_fmpz(L3, lk_callback_z, NULL);
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
  // supply E^2's good primes ascending (incl. p=149 > 142) via Lfunc_use_lpoly_fmpz, then
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
	    uint64_t nmaxO = Lfunc_nmax(LO);
	    // good primes (everything in ap37 except the bad prime 37), ascending,
	    // each raised to k=2; this includes p=149 > nmax(37.a)=142.
	    for (size_t i = 0; i < sizeof(ap37)/sizeof(ap37[0]); i++) {
	      uint64_t p = (uint64_t)ap37[i][0];
	      if (p == 37 || p > nmaxO) continue;
	      eO |= use_power_factor_fmpz(LO, p, 2);
	    }
	    // bad prime 37 supplied LAST (out of order), (1+T)^2:
	    if (37 <= nmaxO)
	      eO |= use_power_factor_fmpz(LO, 37, 2);
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

  // ---- sparse supply: an aborted assembly must be FATAL ----
  // L = M^2 (degree 4, conductor 196 = 14^2, EC-squared mus) supplied with a
  // SINGLE exact square factor at p = 1009, far above nmax(M) for the would-be
  // factor M (degree 2, conductor 14, nmax ~ 100). Certification finds k = 2,
  // but no retained prime is <= nmax(M), so assembly cannot proceed at all.
  // That abort must surface as fatal ERR_POWER: a bare ERR_INSUFF_EULER warning
  // would leave fatal_error() false while rank/zeros/Taylor were never set.
  {
    Lparams_t LpSp = { .degree = 4, .conductor = 196, .normalisation = 0.5,
                       .mus = (double[]){0,0,1,1}, .target_prec = 100, .wprec = 0,
                       .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                       .extract_powers = YES };
    Lerror_t eSp = ERR_SUCCESS;
    Lfunc_t LSp = Lfunc_init_advanced(&LpSp, &eSp);
    assert(!fatal_error(eSp));
    (void) Lfunc_nmax(LSp);
    fmpz_poly_t rSp, fSp;
    fmpz_poly_init(rSp); fmpz_poly_init(fSp);
    fmpz_poly_set_coeff_ui(rSp, 0, 1);           // r = 1 - T + 1009*T^2
    fmpz_poly_set_coeff_si(rSp, 1, -1);
    fmpz_poly_set_coeff_ui(rSp, 2, 1009);
    fmpz_poly_pow(fSp, rSp, 2);                  // F_1009 = r^2, full degree 4
    eSp |= Lfunc_use_lpoly_fmpz(LSp, 1009, fSp);
    fmpz_poly_clear(rSp); fmpz_poly_clear(fSp);
    eSp |= Lfunc_compute(LSp);
    assert(fatal_error(eSp));                    // the abort must be fatal...
    assert(eSp & ERR_POWER);
    assert(eSp & ERR_INSUFF_EULER);              // ...and keep the why
    assert(Lfunc_factors(LSp, NULL, NULL) == 0); // nothing half-assembled
    Lfunc_clear(LSp);
  }

  printf("extract_power_test: sparse-supply abort is fatal OK\n");

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
  eMu |= Lfunc_use_all_lpolys_fmpz(LMu, lk_callback_z, NULL); // INSUFF_EULER warning is fine
  eMu |= Lfunc_compute(LMu);
  assert(eMu & ERR_POWER);            // mus not k-divisible -> refuse to extract
  Lfunc_clear(LMu);

  printf("extract_power_test: I2 non-k-divisible mus OK\n");

  // ---- I3: inexact/acb Euler factors must not be certified (rigor) ----
  // E^2 with one coefficient carrying a tiny error ball. acb supply still drives
  // the guard, but it has no exact fmpz provenance, so extraction must be refused.
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

  // ---- C1: a bad factor that is NOT a k-th power must be rejected (soundness) ----
  // Supply 37.a1's GOOD factors as (E_p)^2 but OVERRIDE the bad prime 37 with the
  // non-square 1 - 3T + T^2. The good-prime-only certificate (pre-fix) k-th-rooted
  // this uncertified bad factor blindly and extracted (n_factors=1); the exact-integer
  // certificate over EVERY supplied factor must reject it. Expect ERR_POWER.
  k_param = 2;
  Lparams_t LpBad = { .degree = 4, .conductor = 37u*37u, .normalisation = 0.5,
                      .mus = (double[]){0,0,1,1}, .target_prec = 100, .wprec = 0,
                      .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                      .extract_powers = YES };
  Lerror_t eBad = ERR_SUCCESS;
  Lfunc_t LBad = Lfunc_init_advanced(&LpBad, &eBad);
  assert(!fatal_error(eBad));
  eBad |= Lfunc_use_all_lpolys_fmpz(LBad, lk_badprime_callback_z, NULL); // INSUFF_EULER warning is fine
  eBad |= Lfunc_compute(LBad);
  assert(eBad & ERR_POWER);            // bad factor not a square -> refuse to extract
  Lfunc_clear(LBad);

  printf("extract_power_test: C1 non-square bad factor OK\n");

  // Complex (non-self-dual) boundary: complex acb factors can still trigger the
  // guard, but the exact fmpz certificate cannot be supplied for them. They must
  // therefore be refused under extract_powers=YES.
  {
    Lparams_t LpC = { .degree=4, .conductor=37u*37u, .normalisation=0.5,
                      .mus=(double[]){0,0,1,1}, .target_prec=100, .wprec=0, .gprec=0,
                      .self_dual=NO, .rank=DK, .cache_dir=".", .extract_powers=YES };
    Lerror_t eC = ERR_SUCCESS;
    Lfunc_t LC = Lfunc_init_advanced(&LpC, &eC);
    assert(!fatal_error(eC));
    eC |= Lfunc_use_all_lpolys(LC, lk_complex_callback, NULL);
    eC |= Lfunc_compute(LC);
    assert(eC & ERR_POWER);
    Lfunc_t *fc = NULL; uint64_t *mc = NULL;
    assert(Lfunc_factors(LC, &fc, &mc) == 0);
    Lfunc_clear(LC);
  }
  printf("extract_power_test: complex acb boundary refused OK\n");

  // ---- E^4 (degree 8, even k=4): multi-k candidate recovery (Fix C) ----
  // L = E^4 has conductor 37^4 (maximal perfect-power exponent E=4) and 2nd moment ~16,
  // so the moment estimate k0~4. The candidate loop considers divisors of E that are
  // multiples of k0 -- here only {4} -- and the rigorous certificate confirms k=4: every
  // (E_p)^4 is an exact 4th power and the mus split into equal blocks of 4. The single-k
  // (round(sqrt(moment))) predecessor could only ever have tried one exponent; this proves
  // k=4 is recovered and a degree-8 4th power is assembled end-to-end.
  k_param = 4;
  Lfunc_t Eref4 = Lfunc_init(2, 37, 0.5, mus, &ec);
  ec = ERR_SUCCESS;
  ec |= Lfunc_use_all_lpolys(Eref4, e_callback, NULL);
  ec |= Lfunc_compute(Eref4);
  assert(!fatal_error(ec));

  Lparams_t Lp4 = { .degree = 8, .conductor = 37ull*37ull*37ull*37ull, .normalisation = 0.5,
                    .mus = (double[]){0,0,0,0,1,1,1,1}, .target_prec = 100, .wprec = 0,
                    .gprec = 0, .self_dual = YES, .rank = DK, .cache_dir = ".",
                    .extract_powers = YES };
  Lerror_t e4 = ERR_SUCCESS;
  Lfunc_t L4 = Lfunc_init_advanced(&Lp4, &e4);
  assert(!fatal_error(e4));
  e4 |= Lfunc_use_all_lpolys_fmpz(L4, lk_callback_z, NULL);
  e4 |= Lfunc_compute(L4);
  assert(!fatal_error(e4));

  Lfunc_t *f4 = NULL; uint64_t *m4 = NULL;
  assert(Lfunc_factors(L4, &f4, &m4) == 1 && m4[0] == 4); // recovered k=4
  assert(Lfunc_rank(f4[0]) == rankE);
  assert(Lfunc_rank(L4) == 4*rankE);
  arb_srcptr zE4 = Lfunc_zeros(Eref4, 0);
  arb_srcptr z4  = Lfunc_zeros(L4, 0);
  assert(arb_overlaps((arb_ptr)(z4+0), (arb_ptr)(zE4+0)));
  assert(arb_overlaps((arb_ptr)(z4+1), (arb_ptr)(zE4+0)));
  assert(arb_overlaps((arb_ptr)(z4+2), (arb_ptr)(zE4+0)));
  assert(arb_overlaps((arb_ptr)(z4+3), (arb_ptr)(zE4+0))); // quadrupled
  assert(arb_overlaps((arb_ptr)(z4+4), (arb_ptr)(zE4+1)));
  acb_t s4; acb_init(s4); acb_pow_ui(s4, Lfunc_sign(Eref4), 4, 100);
  assert(acb_overlaps(s4, Lfunc_sign(L4))); acb_clear(s4); // sign(E)^4
  acb_t vL4, vE4, vE4k; acb_init(vL4); acb_init(vE4); acb_init(vE4k);
  Lerror_t sv4 = Lfunc_special_value(vE4, Eref4, 1.5, 0.0)
               | Lfunc_special_value(vL4, L4, 1.5, 0.0);
  assert(!fatal_error(sv4));
  acb_pow_ui(vE4k, vE4, 4, 100);
  assert(acb_overlaps(vL4, vE4k));
  acb_clear(vL4); acb_clear(vE4); acb_clear(vE4k);
  { arb_t tk; arb_init(tk);
    arb_pow_ui(tk, Lfunc_Taylor(Eref4), 4, 100);
    assert(arb_overlaps((arb_ptr)Lfunc_Taylor(L4), tk));
    arb_clear(tk); }
  Lfunc_clear(L4); Lfunc_clear(Eref4);

  printf("extract_power_test: E^4 multi-k recovery OK (k=4)\n");

  return 0;
}
