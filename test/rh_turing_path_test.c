// End-to-end coverage of the LFUNC_RH_TURING dispatch arm. Selecting Turing on a
// degree-3 object (Sym^2 of 11.a) routes through turing_check_RH, whose deg>=3
// limitation reports ERR_RH_ERROR, where the Buthe default CONFIRMS. That
// contrast proves the dispatch routes to Turing (it FAILS if the arm is dropped
// or rerouted to Buthe). Certified first zero still asserted. Exit 0 = pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

static long g_ainvs[5] = {0, -1, 1, -10, -20};
static long mod_p(long v, long p) { return ((v % p) + p) % p; }
static long ap_good(long p) {
  long a1=g_ainvs[0],a2=g_ainvs[1],a3=g_ainvs[2],a4=g_ainvs[3],a6=g_ainvs[4], c=1;
  for (long x=0;x<p;x++){ long x2=mod_p(x*x,p), r=mod_p(x2*x,p);
    r=mod_p(r+a2*x2,p); r=mod_p(r+a4*x,p); r=mod_p(r+a6,p);
    for (long y=0;y<p;y++) if (mod_p(y*y+a1*x*y+a3*y,p)==r) c++; }
  return p+1-c;
}
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)prec; (void)pm; acb_poly_zero(poly); long P=(long)p;
  if (P==11) { acb_poly_set_coeff_si(poly,0,1); acb_poly_set_coeff_si(poly,1,-1); return; }
  long a=ap_good(P), t=a*a-P;
  acb_poly_set_coeff_si(poly,0,1); acb_poly_set_coeff_si(poly,1,-t);
  acb_poly_set_coeff_si(poly,2,P*t); acb_poly_set_coeff_si(poly,3,-(P*P*P));
}

int main(void) {
  Lerror_t ec = 0; char cd[] = "."; double m3[3] = {0, 1, 0};
  Lparams_t Lp = {0};
  Lp.degree=3; Lp.conductor=121; Lp.normalisation=1.0; Lp.mus=m3;
  Lp.target_prec=DEFAULT_TARGET_PREC; Lp.self_dual=YES; Lp.rank=DK; Lp.cache_dir=cd;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  assert(Lfunc_set_rh_method(L, LFUNC_RH_TURING) == ERR_SUCCESS);
  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L); // runs the LFUNC_RH_TURING arm
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("Turing-selected(deg3) ecode = 0x%lx\n", (unsigned long)ec);
  assert((ec & ERR_RH_ERROR) != 0); // Turing cannot certify deg 3 (Buthe would) => dispatch discriminator
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  arb_t z0; arb_init(z0);
  arb_set_str(z0, "3.89928149477134478", 100);
  arb_add_error_2exp_si(z0, -30); // broaden to span the computed ball
  assert(arb_overlaps(zeros + 0, z0));
  arb_clear(z0);
  Lfunc_clear(L);
  printf("PASS: LFUNC_RH_TURING dispatch routes to Turing (deg-3 ERR_RH_ERROR).\n");
  return 0;
}
