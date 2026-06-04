// Regression test for the configurable output window (bead lfunctions-31c).
#include <flint/acb_poly.h>
#include <flint/arith.h>
#include <flint/fmpz.h>
#include "glfunc.h"
#include "glfunc_internals.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
#include <cstdint>
#include <sys/stat.h>
using std::map; using std::vector;

// ec_37.a1 Euler factors (degree 2, conductor 37), from examples/ec_37.a1.cpp.
static map<int64_t, vector<int64_t>> ef = {
  {2,{1,2,2}},{3,{1,3,3}},{5,{1,2,5}},{7,{1,1,7}},{11,{1,5,11}},{13,{1,2,13}},
  {17,{1,0,17}},{19,{1,0,19}},{23,{1,-2,23}},{29,{1,-6,29}},{31,{1,4,31}},{37,{1,1}},
  {41,{1,9,41}},{43,{1,-2,43}},{47,{1,9,47}},{53,{1,-1,53}},{59,{1,-8,59}},{61,{1,8,61}},
  {67,{1,-8,67}},{71,{1,-9,71}},{73,{1,1,73}},{79,{1,-4,79}},{83,{1,15,83}},{89,{1,-4,89}},
  {97,{1,-4,97}},{101,{1,-3,101}},{103,{1,-18,103}},{107,{1,12,107}},{109,{1,16,109}},
  {113,{1,18,113}},{127,{1,-1,127}},{131,{1,12,131}},{137,{1,6,137}},{139,{1,-4,139}},
  {149,{1,5,149}},{151,{1,-16,151}},{157,{1,-23,157}},{163,{1,18,163}},{167,{1,12,167}},
  {173,{1,-9,173}},{179,{1,-18,179}},{181,{1,-5,181}},{191,{1,4,191}},
};
static void cb(acb_poly_t poly, uint64_t p, int, int64_t, void*) {
  acb_poly_zero(poly);
  auto it = ef.find((int64_t)p);
  if (it != ef.end())
    for (size_t i = 0; i < it->second.size(); ++i)
      acb_poly_set_coeff_si(poly, i, it->second[i]);
}

// Build ec_37.a1 via the advanced API at a given window and cap, in a private
// cache_dir so runs never poison each other. Returns the accumulated ecode.
static Lfunc_t build_37(double max_t, uint64_t max_fft_NN, const char *cache_dir, Lerror_t *ecode) {
  static double mus[2] = {0,1};
  Lparams_t Lp = {}; // zero-init: future Lparams_t fields default safely
  Lp.degree = 2; Lp.conductor = 37; Lp.normalisation = 0.5; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = DK; Lp.rank = DK; Lp.cache_dir = (char*)cache_dir;
  Lp.max_t = max_t; Lp.max_fft_NN = max_fft_NN;
  // The G cache is only written if cache_dir exists: fopen(...,"w") fails on a
  // missing directory and caching is then SILENTLY skipped. Create it (ignore
  // EEXIST) so the cache is actually exercised by these tests.
  mkdir(cache_dir, 0777);
  *ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_compute(L);
  return L;
}

// init-only ec_37.a1 at a given window/cap. Used to confirm that a rejected
// (NaN / huge / sub-floor) window fails LOUD at Lfunc_init_advanced itself,
// fast, before any Euler factors or compute.
static Lfunc_t init_only_37(double max_t, uint64_t max_fft_NN, const char *cache_dir, Lerror_t *ecode) {
  static double mus[2] = {0,1};
  Lparams_t Lp = {};
  Lp.degree = 2; Lp.conductor = 37; Lp.normalisation = 0.5; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = DK; Lp.rank = DK; Lp.cache_dir = (char*)cache_dir;
  Lp.max_t = max_t; Lp.max_fft_NN = max_fft_NN;
  *ecode = ERR_SUCCESS;
  return Lfunc_init_advanced(&Lp, ecode);
}

// Ramanujan tau: degree 2, conductor 1, motivic weight 11 in the algebraic
// normalisation. As in examples/tau.cpp we run it analytically with
// normalisation 5.5, mus {0,1}, so the analytic mus are {5.5,6.5} and
// mu_max = 6.5. The Euler poly at p is 1 - tau(p) x + p^11 x^2, stored as
// coefficients {1, -tau(p), p^11} (tau(p) via FLINT arith_ramanujan_tau).
// build_tau only needs to reach init: the B <= 0.5+mu_max guard fires from the
// conductor/normalisation/mus alone, before any Euler factor is consumed.
static void tau_cb(acb_poly_t poly, uint64_t p, int, int64_t, void*) {
  acb_poly_zero(poly);
  fmpz_t n, t, pk;
  acb_t c;
  fmpz_init(n); fmpz_init(t); fmpz_init(pk); acb_init(c);
  fmpz_set_ui(n, p);
  arith_ramanujan_tau(t, n);          // tau(p)
  fmpz_neg(t, t);                     // -tau(p)
  fmpz_set_ui(pk, p);
  fmpz_pow_ui(pk, pk, 11);            // p^11
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_set_fmpz(c, t);  acb_poly_set_coeff_acb(poly, 1, c);
  acb_set_fmpz(c, pk); acb_poly_set_coeff_acb(poly, 2, c);
  fmpz_clear(n); fmpz_clear(t); fmpz_clear(pk); acb_clear(c);
}

// Build tau via the advanced API at a given window. init may already reject the
// window (the B/mu floor); only if init succeeds do we supply factors. Returns
// the accumulated ecode; on a too-small window it is fatal straight out of init.
static Lfunc_t build_tau(double max_t, const char *cache_dir, Lerror_t *ecode) {
  static double mus[2] = {0,1};
  Lparams_t Lp = {};
  Lp.degree = 2; Lp.conductor = 1; Lp.normalisation = 5.5; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = DK; Lp.rank = DK; Lp.cache_dir = (char*)cache_dir;
  Lp.max_t = max_t; Lp.max_fft_NN = 0;
  mkdir(cache_dir, 0777); // see build_37: cache is skipped if the dir is absent
  *ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_use_all_lpolys(L, tau_cb, NULL);
  return L;
}

int main() {
  // Task 1 assertion: the new fields exist, defaults (sentinels) reproduce the
  // known ec_37.a1 results computed via the plain Lfunc_init path.
  Lerror_t ec;
  Lfunc_t L = build_37(0.0, 0, "build/wt_cache_default", &ec);
  assert(!fatal_error(ec));
  assert(Lfunc_rank(L) == 1);
  arb_t z0ref; arb_init(z0ref);
  arb_set_str(z0ref, "5.0031700140066586953", 300);
  arb_add_error_2exp_si(z0ref, -50);
  assert(arb_overlaps(Lfunc_zeros(L,0) + 0, z0ref));
  arb_clear(z0ref);
  printf("task1 ok\n");

  // Task 2: an explicit max_t equal to the default reproduces the sentinel result.
  Lerror_t ec2;
  Lfunc_t Le = build_37(64.0/2.0, 0, "build/wt_cache_expl", &ec2); // H = 64/degree
  assert(!fatal_error(ec2));
  assert(Lfunc_rank(Le) == 1);
  for (int i = 0; i < 10; ++i)
    assert(arb_overlaps(Lfunc_zeros(Le,0)+i, Lfunc_zeros(L,0)+i));
  Lfunc_clear(Le);
  printf("task2 ok\n");

  // Task 3a: enlarging needs a bigger transform; with the default cap it fails loudly.
  // H=48 requires want_fft_NN=2^17 > default cap 2^16, so ERR_WINDOW_TOO_LARGE fires.
  Lerror_t ec3;
  Lfunc_t Lbig_fail = build_37(48.0, 0 /*cap=2^16*/, "build/wt_cache_toobig", &ec3);
  assert(ec3 & ERR_WINDOW_TOO_LARGE);
  assert(fatal_error(ec3));
  if (Lbig_fail) Lfunc_clear(Lbig_fail);

  // Task 3b: raise the cap to 2^18 and enlarge from H=32 to H=64; the low zeros are
  // unchanged and strictly more zeros are found. H=64 needs fft_NN=2^17, and its G grid
  // exceeds the historical fft_N=2048, so this exercises the fft_N scaling: if fft_N is
  // too small the cyclic convolution aliases and the low zeros shift, breaking the overlap.
  Lerror_t ec3b;
  Lfunc_t Lbig = build_37(64.0, (uint64_t)1<<18, "build/wt_cache_enlarge", &ec3b);
  assert(!fatal_error(ec3b));
  // |eps| = 1
  { arb_t m; arb_init(m); acb_abs(m, Lfunc_sign(Lbig), 100);
    arb_sub_ui(m, m, 1, 100); assert(arb_contains_zero(m)); arb_clear(m); }
  // every default zero (all < 32) reappears in the enlarged run
  int ndef = 0; while (ndef < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(L,0)+ndef)) ndef++;
  int nbig = 0; while (nbig < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(Lbig,0)+nbig)) nbig++;
  assert(nbig > ndef);
  for (int i = 0; i < ndef; ++i)
    assert(arb_overlaps(Lfunc_zeros(Lbig,0)+i, Lfunc_zeros(L,0)+i));
  Lfunc_clear(Lbig);
  printf("task3 ok\n");

  // ---- bead 31c.6: a cache written for one window must not poison another ----
  //
  // The on-disk G-cache filename is keyed only on mus, so two windows with the
  // same mus and default precision collide on ONE file. read_gheader must also
  // validate the window (one_over_B): otherwise a cache written for a LARGER
  // window (B=512) is silently reused for a SMALLER one (B=256) -- the body
  // overwrites L->one_over_B with the cached B and the whole L-function is
  // computed on the wrong grid. The ordering is load-bearing: the larger window
  // has a higher gprec, so writing it FIRST makes the cached gprec pass the
  // (cached >= required) sufficiency check when the smaller window reuses it,
  // unmasking the window bug specifically. On fixed code the one_over_B header
  // mismatch makes the file STALE -> recompute+overwrite at B=256 -> correct.
  {
    const char *shared = "build/wt_cache_poison";
    // 1) Enlarged: H=64 (B=512), cap 2^18. Writes the poisoning cache file.
    Lerror_t ecp;
    Lfunc_t Lp_big = build_37(64.0, (uint64_t)1<<18, shared, &ecp);
    assert(!fatal_error(ecp));
    // The cache must actually have been written (dir exists -> fopen "w" ok).
    {
      struct stat sb;
      assert(stat("build/wt_cache_poison/g_0.5_1.5", &sb) == 0);
    }
    Lfunc_clear(Lp_big);

    // 2) Default window (H=32, B=256) in the SAME directory. On BROKEN code
    // compute_g reads the enlarged file and the body OVERWRITES L->one_over_B to
    // 1/512 (B=512), so the rest of init derives A=fft_NN/512=128 (half the
    // correct 256) and the whole transform runs on the wrong geometry; on FIXED
    // code the one_over_B header mismatch makes the file STALE -> recompute at
    // B=256.
    Lerror_t ecp2;
    Lfunc_t Lp_def = build_37(0.0, 0, shared, &ecp2);
    assert(!fatal_error(ecp2));
    assert(Lfunc_rank(Lp_def) == 1);
    // The low zeros are intrinsic and come out at the right t in both cases, so
    // the first-zero VALUE alone does not distinguish broken from fixed: it is
    // only its radius and the window EXTENT that the wrong B corrupts. We assert
    // the extent. A correct B=256 default trusts t up to the Turing zone
    // 1.5*H = 48, so its largest returned zero is < 50 (here ~47.57). The
    // poisoned B=512 geometry trusts twice as far and returns zeros out to ~79,
    // far above 50. Assert no returned zero exceeds 55: FAILS on broken
    // (max ~79.15), PASSES on fixed (max ~47.57).
    arb_t z0ref; arb_init(z0ref);
    arb_set_str(z0ref, "5.0031700140066586953", 300);
    arb_add_error_2exp_si(z0ref, -50);
    assert(arb_overlaps(Lfunc_zeros(Lp_def,0) + 0, z0ref));
    arb_clear(z0ref);
    arb_t hi; arb_init(hi); arb_set_ui(hi, 55);
    int nzp = 0;
    while (nzp < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(Lp_def,0)+nzp)) {
      // every trusted zero of the correctly-windowed default run is below the
      // 1.5*H = 48 Turing reach, hence < 55. A poisoned B=512 run returns zeros
      // up to ~79, violating this on at least one index.
      assert(arb_lt(Lfunc_zeros(Lp_def,0)+nzp, hi));
      nzp++;
    }
    arb_clear(hi);
    Lfunc_clear(Lp_def);
  }
  printf("poison ok\n");

  // ---- bead 31c.4: every rejected/degenerate window must fail LOUD and fast ----

  // (B3) A NaN max_t must NOT silently fall through to the default window. It is
  // a caller error and must be fatal (not a successful default-window build).
  {
    Lerror_t ecn;
    Lfunc_t Ln = init_only_37(NAN, 0, "build/wt_cache_nan", &ecn);
    assert(fatal_error(ecn));
    if (Ln) Lfunc_clear(Ln);
  }
  printf("nan ok\n");

  // (B3) A huge max_t overflows the required sample count; it must be fatal
  // ERR_WINDOW_TOO_LARGE and must return immediately (no hang / no wraparound).
  {
    Lerror_t ech;
    Lfunc_t Lh = init_only_37(1e17, 0, "build/wt_cache_huge", &ech);
    assert(ech & ERR_WINDOW_TOO_LARGE);
    assert(fatal_error(ech));
    if (Lh) Lfunc_clear(Lh);
  }
  printf("huge ok\n");

  // (Lower-bound 1) A tiny max_t at degree 2 yields want_fft_NN below the fft_N
  // floor (1<<11); that is a heap-overflow geometry and must be fatal
  // ERR_WINDOW_TOO_SMALL. H=0.5 => next_pow2(1024*2*0.5)=1024 < 2048.
  {
    Lerror_t ect;
    Lfunc_t Lt = init_only_37(0.5, 0, "build/wt_cache_tiny", &ect);
    assert(ect & ERR_WINDOW_TOO_SMALL);
    assert(fatal_error(ect));
    if (Lt) Lfunc_clear(Lt);
  }
  printf("tiny ok\n");

  // (Lower-bound 2) A high-motivic-weight object where B <= 0.5 + mu_max: tau has
  // analytic mu_max = 6.5, so B = 8H must exceed 7.0 (H > 0.875). H=0.8 gives
  // B=6.4: it clears the fft_N floor (next_pow2(1024*2*0.8)=2048) but fails the
  // mu/beta floor and must be fatal ERR_WINDOW_TOO_SMALL.
  {
    Lerror_t ecm;
    Lfunc_t Lm = build_tau(0.8, "build/wt_cache_tau", &ecm);
    assert(ecm & ERR_WINDOW_TOO_SMALL);
    assert(fatal_error(ecm));
    if (Lm) Lfunc_clear(Lm);
  }
  // (Lower-bound 2b) The beta-preflight branch, distinct from B <= 0.5+mu_max:
  // tau at H=0.91 gives B=7.28 > 7.0, clearing the B > 0.5+mu_max double guard, but
  // the ftwiddle beta interval is non-positive, so the side-effect-free beta preflight
  // must still reject it (ERR_WINDOW_TOO_SMALL). Locks in spec test (e).
  {
    Lerror_t ecb;
    Lfunc_t Lb = build_tau(0.91, "build/wt_cache_tau_b", &ecb);
    assert(ecb & ERR_WINDOW_TOO_SMALL);
    assert(fatal_error(ecb));
    if (Lb) Lfunc_clear(Lb);
  }
  printf("tau ok\n");

  Lfunc_clear(L);
  return 0;
}
