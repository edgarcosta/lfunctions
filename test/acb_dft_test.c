/*
   Pins the FFT sign/normalisation convention used by src/compute.c after the
   switch to FLINT's acb_dft, and guards the numerical tightness that the
   hand-rolled radix-2 FFT used to provide.

   Convention.  FLINT's *forward* acb_dft has kernel e(-jk/n) (NO 1/n); Part A pins it
   and the cyclic convolution (which carries 1/n exactly once) against golden values at
   n=4.  The project's old "acb_ifft" was the *unnormalised inverse* DFT (a forward FFT
   then index-reverse, kernel e(+jk/n)) -- the opposite-sign transform, not FLINT's
   forward one.  src/compute.c's final iFFT can still use the forward transform because
   the data it feeds is Hermitian (do_pre_iFFT_errors: res[NN-k]=conj(res[k])), and on
   Hermitian input the forward and inverse DFTs coincide and are real-valued.

   Tightness (upstream flintlib/flint#2709).  FLINT's rad2 root table accumulates ~log2(n) bits of error, which
   at fft_NN = 2^16 makes a naive plan ~1000x looser than the old FFT (which
   computed every twiddle directly).  Building the plan at wprec + DFT_PLAN_EXTRA_PREC
   restores tightness (the old hand-rolled FFT was empirically ~0.5-0.6x, i.e. the
   plan-margined FLINT is at least as tight).  Part B asserts the margin delivers
   that improvement and that the resulting balls are tight in absolute terms, so a
   future regression (dropping the margin) fails loudly.
*/

#include <assert.h>
#include <flint/acb.h>
#include <flint/acb_dft.h>
#include "glfunc_internals.h" /* DFT_PLAN_EXTRA_PREC */

static acb_t *vec_alloc(slong n) {
  acb_t *v = (acb_t *)malloc(sizeof(acb_t) * n);
  for (slong i = 0; i < n; i++)
    acb_init(v[i]);
  return v;
}
static void vec_free(acb_t *v, slong n) {
  for (slong i = 0; i < n; i++)
    acb_clear(v[i]);
  free(v);
}
static void vec_set2(acb_t *v, slong i, slong re, slong im) {
  arb_set_si(acb_realref(v[i]), re);
  arb_set_si(acb_imagref(v[i]), im);
}
/* deterministic exact-integer fill: zero input radius, so the output radius is
   purely the transform's own accumulated rounding error */
static void fill_pattern(acb_t *v, slong n) {
  for (slong i = 0; i < n; i++)
    vec_set2(v, i, (i % 19) - 9, (i % 13) - 6);
}
/* largest real/imag ball radius over the vector, as a double */
static double vec_max_rad(acb_t *v, slong n) {
  arb_t m, r;
  arb_init(m);
  arb_init(r);
  for (slong i = 0; i < n; i++) {
    arb_get_rad_arb(r, acb_realref(v[i]));
    if (arb_gt(r, m))
      arb_set(m, r);
    arb_get_rad_arb(r, acb_imagref(v[i]));
    if (arb_gt(r, m))
      arb_set(m, r);
  }
  double d = arf_get_d(arb_midref(m), ARF_RND_NEAR);
  arb_clear(m);
  arb_clear(r);
  return d;
}
static void check_golden(acb_t *got, const slong *re, const slong *im, slong n, const char *what) {
  acb_t want;
  acb_init(want);
  for (slong i = 0; i < n; i++) {
    arb_set_si(acb_realref(want), re[i]);
    arb_set_si(acb_imagref(want), im[i]);
    if (!acb_overlaps(got[i], want)) {
      flint_printf("FAIL %s[%wd]: got ", what, i);
      acb_printd(got[i], 18);
      flint_printf("  want %wd%+wd i\n", re[i], im[i]);
      assert(0);
    }
  }
  acb_clear(want);
}

/* Part A: convention -- forward DFT and cyclic convolution vs golden (n = 4).
   These are exactly the FLINT calls src/compute.c uses. */
static void part_A_golden(slong prec) {
  const slong vre[4] = {1, 2, 3, 4}, vim[4] = {0, 0, 0, 0};
  const slong gre[4] = {5, 6, 7, 8}, gim[4] = {0, 0, 0, 0};
  const slong fwd_re[4] = {10, -2, -2, -2}, fwd_im[4] = {0, 2, 0, -2}; /* F-[1,2,3,4] */
  const slong cnv_re[4] = {66, 68, 66, 60}, cnv_im[4] = {0, 0, 0, 0};  /* cyclic conv */

  acb_dft_rad2_t plan;
  acb_dft_rad2_init(plan, 2, prec + DFT_PLAN_EXTRA_PREC); /* n = 1<<2 = 4 */

  acb_t *v = vec_alloc(4);
  for (slong i = 0; i < 4; i++)
    vec_set2(v, i, vre[i], vim[i]);
  acb_dft_rad2_precomp_inplace((acb_ptr)v, plan, prec);
  check_golden(v, fwd_re, fwd_im, 4, "FLINT forward acb_dft");

  acb_t *f = vec_alloc(4), *g = vec_alloc(4), *w = vec_alloc(4);
  for (slong i = 0; i < 4; i++) {
    vec_set2(f, i, vre[i], vim[i]);
    vec_set2(g, i, gre[i], gim[i]);
  }
  acb_dft_convol_rad2_precomp((acb_ptr)w, (acb_srcptr)f, (acb_srcptr)g, 4, plan, prec);
  check_golden(w, cnv_re, cnv_im, 4, "FLINT acb_dft_convol_rad2");

  vec_free(v, 4);
  vec_free(f, 4);
  vec_free(g, 4);
  vec_free(w, 4);
  acb_dft_rad2_clear(plan);
}

/* Part B: tightness. Compare a plan at wprec (lo) vs wprec+DFT_PLAN_EXTRA_PREC (hi).
   Both enclose the same value; the margined plan must be much tighter, proving the
   margin is necessary and effective. */
static void part_B_tightness(int e, slong prec) {
  slong n = (slong)1 << e;
  acb_dft_rad2_t lo, hi;
  acb_dft_rad2_init(lo, e, prec);
  acb_dft_rad2_init(hi, e, prec + DFT_PLAN_EXTRA_PREC);

  /* transform (in place, so each plan gets its own copy of the input) */
  acb_t *a_lo = vec_alloc(n), *a_hi = vec_alloc(n);
  fill_pattern(a_lo, n);
  fill_pattern(a_hi, n);
  acb_dft_rad2_precomp_inplace((acb_ptr)a_lo, lo, prec);
  acb_dft_rad2_precomp_inplace((acb_ptr)a_hi, hi, prec);
  for (slong i = 0; i < n; i++)
    assert(acb_overlaps(a_lo[i], a_hi[i])); /* same true value enclosed */
  double tr_lo = vec_max_rad(a_lo, n), tr_hi = vec_max_rad(a_hi, n);

  /* cyclic convolution (non-destructive, so inputs are shared) */
  acb_t *cf = vec_alloc(n), *cg = vec_alloc(n), *c_lo = vec_alloc(n), *c_hi = vec_alloc(n);
  fill_pattern(cf, n);
  fill_pattern(cg, n);
  acb_dft_convol_rad2_precomp((acb_ptr)c_lo, (acb_srcptr)cf, (acb_srcptr)cg, n, lo, prec);
  acb_dft_convol_rad2_precomp((acb_ptr)c_hi, (acb_srcptr)cf, (acb_srcptr)cg, n, hi, prec);
  for (slong i = 0; i < n; i++)
    assert(acb_overlaps(c_lo[i], c_hi[i]));
  double cr_lo = vec_max_rad(c_lo, n), cr_hi = vec_max_rad(c_hi, n);

  flint_printf("  n=2^%-2d  transform radius %.3g -> %.3g   convolution radius %.3g -> %.3g\n",
               e, tr_lo, tr_hi, cr_lo, cr_hi);

  /* margined plan is at least as tight */
  assert(tr_hi <= tr_lo);
  assert(cr_hi <= cr_lo);
  /* at the large size the margin is essential: it buys >=50x (empirically ~2000x) */
  if (e >= 16) {
    assert(tr_hi * 50.0 < tr_lo);
    assert(cr_hi * 50.0 < cr_lo);
  }
  /* and the margined balls are tight in absolute terms (catches a broken plan even
     if lo were somehow tight too); generous vs the observed ~1.7e-43 / ~1.7e-40 */
  assert(tr_hi < 1e-40);
  assert(cr_hi < 1e-37);

  vec_free(a_lo, n);
  vec_free(a_hi, n);
  vec_free(cf, n);
  vec_free(cg, n);
  vec_free(c_lo, n);
  vec_free(c_hi, n);
  acb_dft_rad2_clear(lo);
  acb_dft_rad2_clear(hi);
}

int main(void) {
  flint_printf("acb_dft_test: convention + tightness\n");
  part_A_golden(200);
  part_B_tightness(11, 165); /* fft_N  = 2^11 */
  part_B_tightness(16, 165); /* fft_NN = 2^16 */
  flint_printf("acb_dft_test: PASSED\n");
  return 0;
}
