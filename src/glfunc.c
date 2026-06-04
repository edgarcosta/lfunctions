#include "glfunc.h"
#include "assert.h"
#include "glfunc_internals.h"
#include "math.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

static uint64_t next_pow2_u64(uint64_t n) {
  uint64_t p = 1;
  while (p < n) {
    if (p > (UINT64_MAX >> 1)) return p; // never overflow to 0 / loop forever on absurd n
    p <<= 1;
  }
  return p;
}

// Free everything Lfunc_init_advanced allocated before a window guard rejects, then
// signal the error (origin's free-on-error convention; have_B selects whether L->B
// was initialised yet). The caller receives NULL and never double-frees.
static Lfunc_t window_reject(Lfunc *L, bool have_B, Lerror_t *ecode, Lerror_t code) {
  if (have_B) arb_clear(L->B);
  arb_clear(L->pi);
  arb_clear(L->zero_error);
  arb_clear(L->zero_prec);
  free(L->mus);
  free(L);
  *ecode |= code;
  return (Lfunc_t)NULL;
}

bool fatal_error(Lerror_t ecode) { return ecode & 0xFFFFFFFF; }

void fprint_errors(FILE *f, Lerror_t ecode) {
  // fatal errors
  if (ecode & ERR_OOM)
    fprintf(f, "Ran out of memory somewhere.\n");
  if (ecode & ERR_M_ERROR)
    fprintf(f, "Error computing M. Failed Lemma 2 of M_error.tex.\n");
  if (ecode & ERR_NO_DATA)
    fprintf(f, "Looks like we have no usable data.\n");
  if (ecode & ERR_ZERO_ERROR)
    fprintf(f, "Fatal error looking for zeros.\n");
#ifdef BUTHE
  if (ecode & ERR_BUT_ERROR)
    fprintf(
        f,
        "Error doing Buthe check. Estimate for Wf+Winf-Ws* must allow >=0.\n");
#endif
  if (ecode & ERR_UPSAMPLE)
    fprintf(f, "Error computing bounds for upsampling.\n");
  if (ecode & ERR_MU_HALF)
    fprintf(f, "We expect all mu's to be 1/2 non-negative integers.\n");
  if (ecode & ERR_STAT_POINT)
    fprintf(f, "Fatal error in stationary point routine.\n");
  if (ecode & ERR_SPEC_VALUE)
    fprintf(f, "Fatal error in special value routine.\n");
  // warnings
  if (ecode & ERR_INSUFF_EULER)
    fprintf(f, "Don't appear to have enough Euler factors.\n");

  if (ecode & ERR_SOME_DATA)
    fprintf(f, "Data became unusable after output region but before end of "
               "Turing zone.\n");
  if (ecode & ERR_ZERO_PREC)
    fprintf(f, "Couldn't isolate all zeros to requested precision.\n");
  if (ecode & ERR_NO_RANK)
    fprintf(f, "Could not determine rank of L.\n");
  if (ecode & ERR_CONFLICT_RANK)
    fprintf(f, "Computed rank did not agree with what we were told.\n");
  if (ecode & ERR_RH_ERROR)
    fprintf(f, "Failed to confirm RH for zeros in output region.\n");
  if (ecode & ERR_DBL_ZERO)
    fprintf(
        f,
        "Stationary point routine failed to converge. Possible double zero?\n");
  if (ecode & ERR_SPEC_PREC)
    fprintf(
        f, "Failed to achieve desired error bound in Special Value routine.\n");
  if (ecode & ERR_G_INFILE)
    fprintf(f, "Problem opening cached G data file.\n");
  if (ecode & ERR_G_OUTFILE)
    fprintf(f, "Problem opening file to cache G data.\n");
  if (ecode & ERR_G_EXTENT)
    fprintf(f,
            "G data does not extend low enough for this conductor (the cached "
            "or computed grid floor is too high).\n");
  if (ecode & ERR_WINDOW_TOO_LARGE)
    fprintf(f, "Requested output window too large for max_fft_NN (raise max_fft_NN).\n");
  if (ecode & ERR_WINDOW_TOO_SMALL)
    fprintf(f, "Requested output window too small for a valid error analysis.\n");
  if (ecode & ERR_BAD_DEGREE)
    fprintf(f, "The degree of the L-function must be between 1 and %d\n",
            MAX_DEGREE + 1);
  if (ecode & ERR_SPEC_NZ)
    fprintf(f, "Special value routine only works for Im s>=0.\n");
}

// what is the decay (in bits) in the gamma factors from 1/2 to 1/2+i(64/r)
uint64_t decay(Lfunc *L) {
  arb_t tmp1, tmp2, tmp3;
  acb_t s;
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  acb_init(s);
  arb_set_d(acb_realref(s), 0.5);
  abs_gamma(tmp1, s, L, 100);
  arb_set_d(acb_imagref(s), L->max_t);
  abs_gamma(tmp2, s, L, 100);
  arb_div(tmp3, tmp1, tmp2, 100);
  if (verbose) {
    printf("Gamma factor decay is ");
    arb_printd(tmp3, 20);
    printf("\n");
  }
  arb_set_ui(tmp1, 1);
  for (uint64_t j = 0;; j++) {
    if (arb_lt(tmp3, tmp1)) {
      arb_clear(tmp1);
      arb_clear(tmp2);
      arb_clear(tmp3);
      acb_clear(s);
      return j;
    }
    arb_mul_2exp_si(tmp1, tmp1, 1);
  }
}

bool is_half_int(double x) {
  return (x >= 0.0) && ((2.0 * x) == ceil(2.0 * x)) &&
         ((2.0 * x) == floor(2.0 * x));
}

int double_comp(const void *a, const void *b) {
  double *x = (double *)a;
  double *y = (double *)b;
  if (*x < *y)
    return -1;
  else if (*x > *y)
    return 1;
  else
    return 0;
}

Lfunc_t Lfunc_init_advanced(Lparams_t *Lp, Lerror_t *ecode) {
  ecode[0] = ERR_SUCCESS;
  uint64_t i, j;
  arb_t tmp;

  if ((Lp->degree < 2) || (Lp->degree > MAX_DEGREE)) {
    ecode[0] |= ERR_BAD_DEGREE;
    return ((Lfunc_t)NULL);
  }

  Lfunc *L = (Lfunc *)malloc(sizeof(Lfunc));
  if (!L) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  L->degree = Lp->degree;
  L->normalisation = Lp->normalisation;
  L->conductor = Lp->conductor;
  L->mus = (double *)malloc(sizeof(double) * L->degree);
  if (!L->mus) {
    ecode[0] |= ERR_OOM;
    free(L);
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->degree; i++) {
    L->mus[i] = Lp->mus[i] + Lp->normalisation; // alg->anal
    if (!is_half_int(L->mus[i])) {
      ecode[0] |= ERR_MU_HALF;
      // only L and L->mus are live here (no arb_t initialised yet)
      free(L->mus);
      free(L);
      return (Lfunc_t)NULL;
    }
  }

  // we sort the mus so we can name G cache files canonically
  qsort(L->mus, L->degree, sizeof(double), double_comp);

  L->target_prec = Lp->target_prec;
  arb_init(L->zero_prec);
  arb_set_ui(L->zero_prec, 1);
  arb_mul_2exp_si(L->zero_prec, L->zero_prec, -L->target_prec - 1);
  arb_init(L->zero_error);
  arb_add_error(L->zero_error, L->zero_prec);
  arb_init(L->pi);

  // Resolve window geometry from H (max_t) before decay() reads L->max_t.
  // A returned interval that fails to contain the true value is the worst
  // outcome, so every rejected/degenerate window must fail with a fatal
  // Lerror_t here, never silently fall through to a usable-looking result.
  L->max_fft_NN = (Lp->max_fft_NN > 0) ? Lp->max_fft_NN : ((uint64_t)1 << 16);
  uint64_t want_fft_NN;
  // A supplied (non-sentinel) max_t must be finite and strictly positive.
  // A non-finite (NaN/Inf) value is a caller error: it must NOT silently take
  // the default-window branch (the bare Lp->max_t > 0.0 test let NaN through,
  // since NaN > 0.0 is false). The sentinel is exactly 0.0, which is finite.
  if (!isfinite(Lp->max_t)) {
    return window_reject(L, false, ecode, ERR_WINDOW_TOO_SMALL);
  }
  if (Lp->max_t > 0.0) {                       // non-default window
    L->max_t = Lp->max_t;
    L->one_over_B = 1.0 / ((double)OUTPUT_RATIO * L->max_t);
    // Compute the required sample count in binary64 first and reject an
    // overflow before the (uint64_t) cast, so an absurd H can never wrap to a
    // small or zero count.
    double need = ceil(1024.0 * (double)L->degree * L->max_t);
    if (!isfinite(need) || need > (double)UINT64_MAX) {
      return window_reject(L, false, ecode, ERR_WINDOW_TOO_LARGE);
    }
    want_fft_NN = next_pow2_u64((uint64_t)need);
  } else if (Lp->max_t == 0.0) {                // sentinel: default window
    L->max_t = 64.0 / (double)L->degree;
    L->one_over_B = (double)L->degree / 512.0;  // identical to the old g.c:848 value
    want_fft_NN = (uint64_t)1 << 16;            // identical to the old glfunc.c:233 value
  } else {                                      // finite, < 0: too small
    return window_reject(L, false, ecode, ERR_WINDOW_TOO_SMALL);
  }

  // Upper bound: the required transform must fit under the cap.
  if (want_fft_NN > L->max_fft_NN) {
    return window_reject(L, false, ecode, ERR_WINDOW_TOO_LARGE);
  }

  // ---- Lower-bound guards (spec step 4); any failure => fatal, before compute_g ----
  // (1) fft_NN floor: below 1<<11 the error/convolution fill writes res[fft_N-1]
  //     out of the buffer sized fft_NN (heap overflow). fft_N never drops below
  //     1<<11, so want_fft_NN must reach it too.
  if (want_fft_NN < ((uint64_t)1 << 11)) {
    return window_reject(L, false, ecode, ERR_WINDOW_TOO_SMALL);
  }

  // Derive B now (needed by the beta preflight below and reused later for L->B).
  // B = 1/one_over_B = OUTPUT_RATIO * H. L->mus is sorted ascending so the last
  // entry is mu_max.
  double B_eff = 1.0 / L->one_over_B;
  double mu_max = L->mus[L->degree - 1];
  arb_init(L->B);
  {
    arb_t tmpB;
    arb_init(tmpB);
    arb_set_d(tmpB, L->one_over_B);
    arb_inv(L->B, tmpB, 100); // refined to wprec later (after wprec is known)
    arb_clear(tmpB);
  }

  // (3) taylor_terms termination needs B > 4/degree (necessary asymptotic floor;
  //     the hard cap in g.c is the belt-and-suspenders backstop).
  if (!(B_eff > 4.0 / (double)L->degree)) {
    return window_reject(L, true, ecode, ERR_WINDOW_TOO_SMALL);
  }

  // (2) ftwiddle truncation: need B > 0.5 + mu_max AND the resulting decay rate
  //     beta strictly positive, else pre_ftwiddle_error is silently garbage.
  //     Use the shared side-effect-free preflight (same formula as
  //     init_ftwiddle_error) rather than inspecting the damage afterwards.
  arb_const_pi(L->pi, 100); // beta preflight (and decay) need pi
  if (!(B_eff > 0.5 + mu_max) || !ftwiddle_beta_positive(L, 100)) {
    return window_reject(L, true, ecode, ERR_WINDOW_TOO_SMALL);
  }

  if (Lp->wprec > 0)
    L->wprec = Lp->wprec;
  else {
    arb_const_pi(L->pi, 100); // for now, needed by decay()
    // allow enough bits so we will get target_prec at height 1/2+i*max_t
    L->wprec = L->target_prec + decay(L) + EXTRA_BITS;
    if (verbose)
      printf("working precision set to %" PRId64 "\n", L->wprec);
  }
  arb_const_pi(L->pi, L->wprec); // set it properly now we know what wprec is
  L->gprec = Lp->gprec;
  L->self_dual = Lp->self_dual;
  L->rank = Lp->rank;
  L->cache_dir = Lp->cache_dir;

  // See Lemma 2 of M_error1.pdf, Lemma 5 of g.pdf
  // r is always >=2
  L->nus = (arb_t *)malloc(sizeof(arb_t) * L->degree);
  if (!L->nus) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (j = 0; j < 2; j++) {
    arb_init(L->nus[j]);
    arb_set_d(L->nus[j], 0.5 * (L->mus[j] - 0.5));
  }
  for (j = 2; j < L->degree; j++) {
    arb_init(L->nus[j]);
    arb_set_d(L->nus[j], 0.5 * (L->mus[j] - 1.0));
  }

  // See Lemma 2/5 again
  arb_init(L->nu);
  arb_set(L->nu, L->nus[0]);
  for (j = 1; j < L->degree; j++)
    arb_add(L->nu, L->nu, L->nus[j], L->target_prec);
  arb_div_ui(L->nu, L->nu, L->degree, L->target_prec);
  arb_mul_2exp_si(L->nu, L->nu, 1);
  arb_init(tmp);
  arb_set_d(tmp, 0.5);
  arb_add(L->nu, L->nu, tmp, L->target_prec);

  // mu=-1/2+1/r(1+sum mu_j) See Lemma 5.2 of Artin
  arb_init(L->mu);
  arb_set_d(L->mu, -0.5);
  double smu = 1.0;
  for (i = 0; i < L->degree; i++)
    smu += L->mus[i];
  arb_set_d(tmp, smu);
  arb_div_ui(tmp, tmp, L->degree, L->target_prec);
  arb_add(L->mu, L->mu, tmp, L->target_prec);

  ecode[0] |= compute_g(L);
  if (fatal_error(ecode[0]))
    return (Lfunc_t)NULL;
  if (verbose) {
    printf("eq 5-9 error = ");
    arb_printd(L->eq59, 20);
    printf("\n");
  }

  // L->B was already init'd and set at 100 bits for the window preflight above;
  // refine it to full working precision now that wprec is known.
  arb_set_d(tmp, L->one_over_B);
  arb_inv(L->B, tmp, L->wprec);

  arb_init(L->two_pi_by_B);
  arb_set_d(L->two_pi_by_B, L->one_over_B * 2.0);
  arb_mul(L->two_pi_by_B, L->two_pi_by_B, L->pi, L->wprec);

  // fft_N (short Euler-convolution length) must exceed the linear-convolution support
  // of the G grid [low_i,hi_i] and the coefficient buckets, else the cyclic load
  // n%fft_N (compute.c) silently aliases. The lowest bucket is for n=1 at
  // calc_m(1) = floor(log(1/sqrt(N))/(2*pi/B)), which is conductor-DEPENDENT (it
  // sinks as N grows); conv_support() uses that floor (see glfunc_internals.h).
  // max(1<<11,...) keeps every default window bit-for-bit (support <= ~1350 < 2048).
  {
    uint64_t support = conv_support(L); // G grid + coeff buckets (see glfunc_internals.h)
    uint64_t floor_n = (uint64_t)1 << 11;
    L->fft_N = next_pow2_u64(support > floor_n ? support : floor_n);
  }
  L->fft_NN = want_fft_NN; // final output length (= 1<<16 for the default window)
  // Defensive backstop only: this branch is UNREACHABLE for any accepted input.
  // want_fft_NN >= 1<<11 (the sub-floor guard above rejects anything smaller), and
  // fft_N = next_pow2(max(1<<11, support)). support ~ (hi_i-low_i) + (hi_i-calc_m(1))
  // grows like B*(1+log N) = 8*max_t*(1+log N), while fft_NN ~ 1024*degree*max_t; for
  // fft_N to exceed fft_NN we would need 1+log N > ~128*degree, i.e. log N > 255 at
  // degree 2, far beyond any uint64_t conductor (log N <= ~44). So fft_N <= fft_NN
  // always. L is only PARTIALLY constructed here (compute_g ran, but the w/ww/res/
  // zeros/upsampling allocations below have not), so Lfunc_clear would dereference
  // not-yet-allocated pointers; we keep the graceful fatal return instead. Freeing the
  // partial Lfunc for such early init failures is tracked by bead lfunctions-1eb.
  if (L->fft_N > L->fft_NN) { ecode[0] |= ERR_WINDOW_TOO_LARGE; return (Lfunc_t)NULL; }

  L->A = L->fft_NN * L->one_over_B;
  arb_init(L->arb_A);
  arb_set_d(L->arb_A, L->A);
  arf_init(L->arf_A);
  arf_set_d(L->arf_A, L->A);

  arb_init(L->one_over_A);
  arb_inv(L->one_over_A, L->arb_A, L->wprec);
  arf_init(L->arf_one_over_A);
  arf_ui_div(L->arf_one_over_A, 1, L->arf_A, L->wprec, ARF_RND_NEAR);

  L->G = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
  if (!L->G) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->fft_N; i++)
    acb_init(L->G[i]);

  L->eta = 0.0;
  arb_init(L->delta);
  arb_mul_2exp_si(L->delta, L->pi, -1); // pi/2
  arb_set_d(tmp, 1.0 - L->eta);
  arb_mul(L->delta, L->delta, tmp, L->wprec); // (1-eta)pi/2
  arb_init(L->exp_delta);
  arb_neg(tmp, L->delta);
  arb_exp(L->exp_delta, tmp, L->wprec);

  L->w =
      (acb_t *)malloc(sizeof(acb_t) * L->fft_N / 2); // twiddles for little FFT
  if (!L->w) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (i = 0; i < L->fft_N / 2; i++)
    acb_init(L->w[i]);
  acb_initfft(L->w, L->fft_N, L->wprec); // set twiddles for little FFT

  L->ww =
      (acb_t *)malloc(sizeof(acb_t) * L->fft_NN / 2); // twiddles for big FFT
  if (!L->ww) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->fft_NN / 2; i++)
    acb_init(L->ww[i]);
  acb_initfft(L->ww, L->fft_NN, L->wprec); // set twiddles for big FFT

  // space for the zeros once we isolate them
  L->zeros[0] = (arb_t *)malloc(sizeof(arb_t) * MAX_ZEROS);
  L->zeros[1] = (arb_t *)malloc(sizeof(arb_t) * MAX_ZEROS);
  if ((!L->zeros[0]) || (!L->zeros[1])) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < MAX_ZEROS; i++) {
    arb_init(L->zeros[0][i]);
    arb_init(L->zeros[1][i]);
  }

  L->kres = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
  if (!L->kres) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  L->skm = (acb_t **)malloc(sizeof(acb_t *) * L->max_K);
  if (!L->skm) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  uint64_t k, n;
  for (k = 0; k < L->max_K; k++) {
    L->skm[k] = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
    if (!L->skm[k]) {
      arb_clear(tmp);
      ecode[0] |= ERR_OOM;
      return (Lfunc_t)NULL;
    }

    for (n = 0; n < L->fft_N; n++)
      acb_init(L->skm[k][n]);
  }

  L->res = (acb_t *)malloc(sizeof(acb_t) * L->fft_NN);
  if (!L->res) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (n = 0; n < L->fft_N; n++)
    acb_init(L->kres[n]);
  for (n = 0; n < L->fft_NN; n++)
    acb_init(L->res[n]);

  arb_init(L->pre_ftwiddle_error);
  arb_init(L->ftwiddle_error);
  init_ftwiddle_error(L, L->wprec);

  arb_init(L->one_over_root_N);
  arb_init(L->sum_ans);
  acb_init(L->sign);
  acb_init(L->sqrt_sign);
  L->allocated_M = 8192;
  L->ans = (acb_t *)malloc(sizeof(acb_t) * L->allocated_M);
  if (!L->ans) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (size_t i = 0; i < L->allocated_M; ++i)
    acb_init(L->ans[i]);

  arb_clear(tmp);

#ifdef BUTHE
  init_buthe(L, L->wprec); // setup stuff for Buthe zero check
#endif
#ifdef TURING
  arb_init(L->imint);
  arb_init(L->X);
#endif
  L->nmax_called = false; // noone has called nmax yet

  arb_init(L->Lam_d);
  arb_init(L->L_d);

  ecode[0] |= init_upsampling(L);
  if (fatal_error(ecode[0]))
    return (Lfunc_t)NULL;

  // (4) Upsample period-fit (spec step 4, last bullet): the output + Turing +
  // upsampling-guard samples must fit in one output period of fft_NN. This is
  // exactly L->u_no_values <= fft_NN, where init_upsampling has just set
  // u_no_values = fft_NN/OUTPUT_RATIO + fft_NN/TURING_RATIO + 4*u_N*u_stride + 1
  // using the same parameter search as the runtime; reading it here (rather than
  // replicating the search) keeps the two from diverging. Runs before zero-
  // finding (which happens in Lfunc_compute), so no degraded result is returned.
  if (L->u_no_values > L->fft_NN) {
    // L is fully constructed here (init_upsampling succeeded, every u_* field and
    // array is set), so the authoritative teardown applies before we hand back NULL.
    Lfunc_clear((Lfunc_t)L);
    ecode[0] |= ERR_WINDOW_TOO_SMALL;
    return (Lfunc_t)NULL;
  }

  return (Lfunc_t)L;
}

Lfunc_t Lfunc_init(uint64_t degree, uint64_t conductor, double normalisation,
                   const double *mus, Lerror_t *ecode) {
  Lparams_t Lp;
  Lp.degree = degree;
  Lp.conductor = conductor;
  Lp.normalisation = normalisation;
  Lp.mus = (double *)malloc(sizeof(double) * degree);
  if (!Lp.mus) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (size_t i = 0; i < degree; ++i)
    Lp.mus[i] = mus[i];
  Lp.target_prec = DEFAULT_TARGET_PREC;
  Lp.rank = DK;
  Lp.self_dual = DK;
  Lp.cache_dir = ".";
  Lp.gprec = 0; // We will try to do something sensible
  Lp.wprec = 0; // ditto
  Lp.max_t = 0.0;     // sentinel => 64/degree
  Lp.max_fft_NN = 0;  // sentinel => 1<<16

  // Lfunc_init_advanced copies Lp.mus into L->mus, so this scratch array is ours
  // to free; otherwise it leaks (Lp is a stack local that goes out of scope).
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  free(Lp.mus);
  return L;
}

int64_t Lfunc_wprec(Lfunc_t Lf) {
  Lfunc *L;
  L = (Lfunc *)Lf;
  return L->wprec;
}

#ifdef __cplusplus
}
#endif
