/*
  Regression test for the power / repeated-factor guard (ERR_POWER); bead lfunctions-w70.1.
  Asserts on the error bitfield, never on printed output.

    1. L(chi5)^2          (degree 2, cond 25):  genuine square     -> ERR_POWER by default.
    2. same, opted in     (extract_powers):  confirmed a power, M is degree 1 -> ERR_BAD_DEGREE.
    3. L(chi5)*L(chi7)     (degree 2, cond 35):  imprimitive but squarefree -> NOT rejected.
    4. L(chi3)^4          (degree 4, cond 81):  4th power          -> ERR_POWER by default.
*/
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <flint/acb_poly.h>
#include "glfunc.h"

// real (quadratic) Dirichlet characters, value in {-1,0,1}
static int chi3(uint64_t p){ uint64_t r=p%3; if(r==0) return 0; return (r==1)?1:-1; }
static int chi5(uint64_t p){ uint64_t r=p%5; if(r==1||r==4) return 1; if(r==2||r==3) return -1; return 0; }
static int chi7(uint64_t p){ uint64_t r=p%7; if(r==0) return 0; return (r==1||r==2||r==4)?1:-1; }

// (1 - chi(p) T) as an acb_poly
static void linear_factor(acb_poly_t f, int chi){
  acb_poly_one(f);
  if(chi!=0) acb_poly_set_coeff_si(f,1,-chi);
}

// kind: 1 = (1-chi5 T)^2, 2 = (1-chi5 T)(1-chi7 T), 3 = (1-chi3 T)^4
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param){
  (void)d;
  int kind = *(int *)param;
  acb_poly_t a,b;
  acb_poly_init(a); acb_poly_init(b);
  if(kind==1){
    linear_factor(a, chi5(p));
    acb_poly_mul(poly, a, a, prec);                 // (1-chi5 T)^2
  } else if(kind==2){
    linear_factor(a, chi5(p));
    linear_factor(b, chi7(p));
    acb_poly_mul(poly, a, b, prec);                 // (1-chi5 T)(1-chi7 T)
  } else if(kind==3) {
    linear_factor(a, chi3(p));
    acb_poly_mul(b, a, a, prec);                    // (1-chi3 T)^2
    acb_poly_mul(poly, b, b, prec);                 // (1-chi3 T)^4
  } else if(kind==4){
    // synthetic: (1 - 0.5 T)^2 -- a non-squarefree factor with |a_p| = 1, so the 2nd
    // moment stays ~1 (below POWER_MOMENT_THRESHOLD). It is NOT a real arithmetic Euler
    // factor; it exists only to drive the guard's moment signal below threshold so that,
    // at conductor 25 (=5^2), ONLY the conductor perfect-power signal can trigger the reject.
    acb_t half;
    acb_init(half);
    acb_set_d(half, -0.5);
    acb_poly_one(a);
    acb_poly_set_coeff_acb(a, 1, half);   // 1 - 0.5 T
    acb_poly_mul(poly, a, a, prec);        // (1 - 0.5 T)^2
    acb_clear(half);
  }
  acb_poly_clear(a); acb_poly_clear(b);
}

// supply Euler factors and run the computation; return the accumulated error code
static Lerror_t run(Lfunc_t L, int kind){
  Lerror_t ecode = ERR_SUCCESS;
  int k = kind;
  ecode |= Lfunc_use_all_lpolys(L, cb, &k);
  ecode |= Lfunc_compute(L);
  return ecode;
}

int main(void){
  double mus2[] = {0,0};       // L(chi5)^2 -> [0,0]
  double mus_57[] = {0,1};     // L(chi5) even (0), L(chi7) odd (1)
  double mus4[] = {1,1,1,1};   // chi3 is odd (mu=1); chi3^4 -> [1,1,1,1]
  Lerror_t ecode;

  // 1. genuine square -> rejected
  ecode = ERR_SUCCESS;
  Lfunc_t L1 = Lfunc_init(2, 25, 0.0, mus2, &ecode);
  assert(!fatal_error(ecode));
  Lerror_t e1 = run(L1, 1);
  Lfunc_clear(L1);
  assert(e1 & ERR_POWER);

  // 2. same object, opted in: extraction attempted, M would be degree 1 -> ERR_BAD_DEGREE.
  //    The power IS recognised (not ERR_POWER), but the primitive factor is sub-degree-2.
  Lparams_t Lp;
  Lp.degree = 2; Lp.conductor = 25; Lp.normalisation = 0.0; Lp.mus = mus2;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = DK; Lp.rank = DK; Lp.cache_dir = ".";
  Lp.extract_powers = YES;
  ecode = ERR_SUCCESS;
  Lfunc_t L2 = Lfunc_init_advanced(&Lp, &ecode);
  assert(!fatal_error(ecode));
  Lerror_t e2 = run(L2, 1);
  Lfunc_clear(L2);
  assert(e2 & ERR_BAD_DEGREE);   // confirmed a power; M = degree 1 is unsupported
  assert(!(e2 & ERR_POWER));     // not guard-rejected; extraction was attempted

  // 3. imprimitive but squarefree (two distinct primitives) -> NOT rejected
  ecode = ERR_SUCCESS;
  Lfunc_t L3 = Lfunc_init(2, 35, 0.0, mus_57, &ecode);
  assert(!fatal_error(ecode));
  Lerror_t e3 = run(L3, 2);
  Lfunc_clear(L3);
  assert(!(e3 & ERR_POWER));

  // 4. a 4th power at degree 4 -> rejected
  ecode = ERR_SUCCESS;
  Lfunc_t L4 = Lfunc_init(4, 81, 0.0, mus4, &ecode);
  assert(!fatal_error(ecode));
  Lerror_t e4 = run(L4, 3);
  Lfunc_clear(L4);
  assert(e4 & ERR_POWER);

  // 5. synthetic pure square: 2nd moment ~1 (below threshold) but conductor 25 = 5^2.
  //    No squarefree full-degree factor + moment < 3.5, so ONLY the conductor signal
  //    can reject. Fails without the conductor trigger (the moment alone lets it through).
  ecode = ERR_SUCCESS;
  Lfunc_t L5 = Lfunc_init(2, 25, 0.0, mus2, &ecode);
  assert(!fatal_error(ecode));
  Lerror_t e5 = run(L5, 4);
  Lfunc_clear(L5);
  assert(e5 & ERR_POWER);

  printf("power_guard_test: all assertions passed\n");
  return 0;
}
