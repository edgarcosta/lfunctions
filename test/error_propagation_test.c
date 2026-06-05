/*
 * Regression test for the error-propagation bug in src/error.c
 * do_pre_iFFT_errors (bead lfunctions-fxm).
 *
 * Bug: the loop at lines ~537-543 iterated j=0..i-1 over the populated output
 * bins but added the F_hat_twiddle aliasing error and the eq-59 G-truncation
 * error to res[i] (the post-break bin) i times instead of to res[j].  The
 * populated bins res[0..i-1] received neither error term; res[i] received both
 * i times, then was immediately zeroed by acb_zero() in the subsequent loop
 * (line ~623).  The certified radius of every output ball was therefore
 * under-reported.
 *
 * Fix: change the four L->res[i] to L->res[j] in that loop.
 *
 * Why public-API ball radii do NOT distinguish fixed from broken:
 *   fhattwiddle and err (= L->eq59 * L->sum_ans) are ~2^{-wprec} and, on every
 *   reachable configuration, are dominated by M_error on the SAME populated bins.
 *   Lowering gprec inflates eq59 but inflates M_error in lockstep (the G grid
 *   shrinks, so the coefficient cutoff M drops), so the output radius moves
 *   together for both versions -- e.g. at gprec 50 the Taylor radius is ~4.7e-9
 *   either way.  No radius threshold on a public output separates them.
 *
 * The test therefore has two parts:
 *   (1) Correctness/regression (full public compute): rank, first zero, and the
 *       Taylor coefficient overlapping the known BSD constant.
 *   (2) A white-box discriminator that FAILS on the broken res[i] code.  M_error
 *       does not depend on sum_ans, so inflating sum_ans on a fresh handle and
 *       calling do_pre_iFFT_errors directly makes err dominate res[0]: the fix
 *       (res[j]) leaves res[0] with radius ~= eq59*sum_ans (~1e-33), while the bug
 *       (res[i]) sends err to res[break] -- which src/error.c ~683 then zeroes --
 *       leaving res[0] at ~= M_error[0] (~1e-61).
 */

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include <flint/mag.h>
#include "glfunc.h"
#include "glfunc_internals.h"  /* white-box: Lfunc, L->res/sum_ans, do_pre_iFFT_errors */

/* Euler factors for ec_37.a1 (from the LMFDB sidebar).
 * https://www.lmfdb.org/EllipticCurve/Q/37/a/1 */
static const int64_t ef_p[]  = {
  2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37,
  41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139
};
/* a1 coefficient for each prime (a0 = 1 always) */
static const int64_t ef_a1[] = {
   2,  3,  2,  1,  5,  2,  0,  0, -2, -6,  4,  1,
   9, -2,  9, -1, -8,  8, -8, -9,  1, -4, 15, -4, -4,
  -3, -18, 12, 16, 18, -1, 12,  6, -4
};
/* a2 = p for good primes; p=37 is bad (conductor factor) => a2=0 */
static const int64_t ef_a2[] = {
   2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31,  0,
  41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
  101, 103, 107, 109, 113, 127, 131, 137, 139
};
#define N_PRIMES ((int)(sizeof(ef_p)/sizeof(ef_p[0])))

static void lpoly_callback(acb_poly_t poly, uint64_t p,
                           int d __attribute__((unused)),
                           int64_t prec __attribute__((unused)),
                           void *param __attribute__((unused)))
{
  for (int k = 0; k < N_PRIMES; k++) {
    if ((uint64_t)ef_p[k] == p) {
      acb_poly_zero(poly);
      acb_poly_set_coeff_si(poly, 0, 1);
      acb_poly_set_coeff_si(poly, 1, ef_a1[k]);
      if (ef_a2[k] != 0)
        acb_poly_set_coeff_si(poly, 2, ef_a2[k]);
      return;
    }
  }
  /* prime not in table: signal end of factors */
  acb_poly_zero(poly);
}

int main(void)
{
  Lerror_t ecode = ERR_SUCCESS;
  double mus[2] = {0.0, 1.0};

  Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ecode);
  if (fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    return 1;
  }

  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if (fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    Lfunc_clear(L);
    return 1;
  }

  ecode |= Lfunc_compute(L);
  if (fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    Lfunc_clear(L);
    return 1;
  }

  /* Report the Taylor radius (informational). */
  arb_srcptr taylor = Lfunc_Taylor(L);
  double r_d = mag_get_d(arb_radref(taylor));
  printf("Taylor = ");
  arb_printd(taylor, 20);
  printf("\nTaylor radius = %.6e\n", r_d);

  /*
   * Verify the Taylor coefficient overlaps the known BSD constant.
   * This catches any regression introduced by the fix (the certified ball
   * must still contain the truth).
   *
   * NOTE: a radius threshold on this PUBLIC output cannot distinguish the broken
   * res[i] from the fixed res[j] (M_error dominates the omitted terms on the same
   * bins, in lockstep across gprec).  The distinguishing check is the white-box
   * res[0] assertion at the end of this file.
   */
  arb_t bsd;
  arb_init(bsd);
  arb_set_str(bsd,
    "0.305999773834052301820483683321676474452637774590771998534541832481"
    "016050469290169911495257337795897237898682879524967997997869651621709"
    "648704953228700", 400);
  assert(arb_overlaps(taylor, bsd));
  printf("Taylor overlaps BSD constant: PASS\n");
  arb_clear(bsd);

  /* Verify rank and first zero. */
  assert(Lfunc_rank(L) == 1);
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  arb_t ref;
  arb_init(ref);
  arb_set_str(ref, "5.0031700140066586953", 300);
  arb_add_error_2exp_si(ref, -50);
  assert(arb_overlaps(zeros, ref));
  printf("First zero overlaps reference: PASS\n");
  arb_clear(ref);

  /*
   * White-box discriminator (FAILS on the broken res[i] code; see the header).
   * M_error ignores sum_ans, so inflate sum_ans on a fresh handle and call
   * do_pre_iFFT_errors directly: the fix routes err = eq59*sum_ans onto the
   * populated bin res[0]; the bug routes it to res[break], later zeroed.
   */
  {
    Lerror_t e = ERR_SUCCESS;
    Lfunc_t Lw = Lfunc_init(2, 37, 0.5, mus, &e);
    e |= Lfunc_use_all_lpolys(Lw, lpoly_callback, NULL);
    assert(!fatal_error(e));
    Lfunc *Li = (Lfunc *)Lw;
    arb_set_d(Li->sum_ans, 1e30);   /* err = eq59*sum_ans ~ 1e-33, far above any M_error */
    acb_set_ui(Li->res[0], 1);      /* unit midpoint: epsilon branch well-defined, |sqrt_sign|=1 */
    e = do_pre_iFFT_errors(Li);
    assert(!fatal_error(e));
    double r0 = mag_get_d(arb_radref(acb_realref(Li->res[0])));
    printf("res[0] radius after inflated do_pre_iFFT_errors = %.4e\n", r0);
    /* FIXED: err landed on res[0] -> r0 ~ 1e-33.  BROKEN: err went to res[break]
     * and was zeroed -> r0 ~ M_error[0] ~ 1e-61.  Threshold 1e-40 has ~7 orders
     * of margin above the fixed value and ~21 below the broken value. */
    assert(r0 > 1e-40);
    Lfunc_clear(Lw);
    printf("res[j] error-distribution discriminator: PASS\n");
  }

  Lfunc_clear(L);
  printf("PASS: error_propagation_test\n");
  return 0;
}
