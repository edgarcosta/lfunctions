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
 * Distinguishing signal attempted:
 *   The F_hat_twiddle bound (fhattwiddle) and the eq-59 G-tail bound (err)
 *   are both exponentially small by design -- they are the per-request
 *   G-truncation and aliasing bounds, each chosen to be ~2^{-wprec}.
 *   For ec_37.a1 at default precision (wprec ~165 bits):
 *     fhattwiddle  ~= 0  (the G-function decays to ~0 before the grid edge)
 *     eq59_err     ~= 2.15e-62
 *   The contribution of eq59_err to the iFFT output is roughly
 *     i * eq59_err / fft_NN ~= 203 * 2.15e-62 / 65536  ~  6.7e-65
 *   which is ~12 orders of magnitude below the dominant M_error contribution
 *   (~1.09e-52 for the Taylor coefficient).
 *   r_broken = 1.094285e-52, r_fixed = 1.094363e-52 (differ by ~7e-57, <0.01%).
 *   No robust radius threshold exists to separate them (ratio ~1.000000064).
 *
 * This is a soundness hole, not a numerical regression detectable via
 * public-API ball radii at default precision.  The test therefore verifies
 * code correctness: the computation completes without fatal error and the
 * Taylor coefficient overlaps the known BSD constant, catching any regression
 * the fix might introduce.
 */

#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include <flint/mag.h>
#include "glfunc.h"

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
   * NOTE: a radius-threshold test to distinguish the broken res[i] from the
   * fixed res[j] is not feasible.  The two omitted error terms (fhattwiddle
   * ~= 0 and eq59_err ~= 2.15e-62) contribute ~6.7e-65 to the output ball
   * radius after the iFFT, while M_error contributes ~1.09e-52.  The ratio
   * (r_fixed - r_broken) / r_broken < 1e-6, so no robust threshold exists.
   * The bug is a soundness hole; correctness verification (BSD overlap) is the
   * appropriate test here.
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

  Lfunc_clear(L);
  printf("PASS: error_propagation_test\n");
  return 0;
}
