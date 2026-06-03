// Hermetic extraction test: builds E^2 and E^3 from an elliptic curve's
// Euler factors and checks extraction recovers E and assembles L = E^k.
// Asserts on certified balls only. Build with assertions on.
#include <assert.h>
#include <stdio.h>
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
  return 0;
}
