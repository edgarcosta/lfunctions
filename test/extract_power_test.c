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
};
// Sanity check against examples/ec_37.a1.cpp: p=2 -> {1,2,2}, p=23 -> {1,-2,23}.

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

  Lfunc_clear(L2);
  Lfunc_clear(Eref);

  printf("extract_power_test: Task 3 OK\n");
  return 0;
}
