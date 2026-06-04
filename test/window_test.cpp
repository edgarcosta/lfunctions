// Regression test for the configurable output window (bead lfunctions-31c).
#include <flint/acb_poly.h>
#include <flint/arith.h>
#include <flint/fmpz.h>
#include <flint/ulong_extras.h>
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

// Same ec_37.a1, but supply the ANALYTIC Euler-factor coefficients a_m * p^{-m/2}
// instead of the algebraic a_m. Combined with normalisation = 0 and analytic mus
// {0.5,1.5}, this reconstructs the SAME analytic L-function that (algebraic coeffs,
// normalisation = 0.5, mus {0,1}) produces: the library multiplies coefficient m by
// p^{-m*normalisation} (coeff.c), so {a_m, norm 0.5} and {a_m*p^{-m/2}, norm 0} feed
// identical Dirichlet coefficients, and L->mus = mus + norm is {0.5,1.5} either way.
// Used by the equivalent-(mus,normalisation) window invariant below.
static void cb_anal(acb_poly_t poly, uint64_t p, int, int64_t prec, void*) {
  acb_poly_zero(poly);
  auto it = ef.find((int64_t)p);
  if (it == ef.end()) return;
  arb_t lp, sc, tm, half; arb_init(lp); arb_init(sc); arb_init(tm); arb_init(half);
  arb_log_ui(lp, p, prec);
  acb_t c; acb_init(c);
  for (size_t m = 0; m < it->second.size(); ++m) {
    arb_set_d(half, -0.5 * (double)m);
    arb_mul(tm, lp, half, prec);             // -(m/2) log p
    arb_exp(sc, tm, prec);                    // p^{-m/2}
    arb_mul_si(sc, sc, it->second[m], prec);  // a_m * p^{-m/2}
    acb_set_arb(c, sc);
    acb_poly_set_coeff_acb(poly, (slong)m, c);
  }
  arb_clear(lp); arb_clear(sc); arb_clear(tm); arb_clear(half); acb_clear(c);
}

// Like build_37 but with caller-chosen normalisation, mus, self_dual, rank, and
// Euler-factor callback, so the invariants below can vary exactly one axis at a
// fixed window. Returns the accumulated ecode (see build_37 for the cache note).
typedef void (*wt_cbfn)(acb_poly_t, uint64_t, int, int64_t, void*);
static Lfunc_t build_37_full(double max_t, uint64_t max_fft_NN, double normalisation,
                             double *mus, int self_dual, int rank, wt_cbfn f,
                             const char *cache_dir, Lerror_t *ecode) {
  Lparams_t Lp = {};
  Lp.degree = 2; Lp.conductor = 37; Lp.normalisation = normalisation; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = self_dual; Lp.rank = rank; Lp.cache_dir = (char*)cache_dir;
  Lp.max_t = max_t; Lp.max_fft_NN = max_fft_NN;
  mkdir(cache_dir, 0777);
  *ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_use_all_lpolys(L, f, NULL);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_compute(L);
  return L;
}

// Count the (contiguous, zero-terminated) zeros stored on a side.
static int wt_count_zeros(Lfunc_t L, uint64_t side) {
  int n = 0;
  while (n < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(L, side) + n)) n++;
  return n;
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

// Product of the quadratic characters mod 227 and mod 229: a degree-2, self-dual
// L-function of conductor 227*229 = 51983, so M0 = ceil(sqrt(N)/100) = 3 > 1
// (unlike conductor 37, where M0 = 1). chi_227 is odd (227 = 3 mod 4) and chi_229
// even, giving analytic mus {0,1} at normalisation 0. The good-prime Euler factor
// is (1 - chi_227(p) x)(1 - chi_229(p) x); a ramified chi_p(p) = 0 drops its factor.
static void dir_cb(acb_poly_t poly, uint64_t p, int, int64_t, void*) {
  acb_poly_zero(poly);
  long a = n_jacobi((slong)(p % 227), 227);   // chi_227(p) in {-1,0,1}
  long b = n_jacobi((slong)(p % 229), 229);   // chi_229(p)
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_poly_set_coeff_si(poly, 1, -(a + b));
  if (a != 0 && b != 0) acb_poly_set_coeff_si(poly, 2, a * b);
}
static Lfunc_t build_dir(double max_t, uint64_t max_fft_NN, const char *cache_dir, Lerror_t *ecode) {
  static double mus[2] = {0,1};
  Lparams_t Lp = {};
  Lp.degree = 2; Lp.conductor = 227 * 229; Lp.normalisation = 0.0; Lp.mus = mus;
  Lp.target_prec = DEFAULT_TARGET_PREC; Lp.wprec = 0; Lp.gprec = 0;
  Lp.self_dual = DK; Lp.rank = DK; Lp.cache_dir = (char*)cache_dir;
  Lp.max_t = max_t; Lp.max_fft_NN = max_fft_NN;
  mkdir(cache_dir, 0777);
  *ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_use_all_lpolys(L, dir_cb, NULL);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_compute(L);
  return L;
}

int main() {
  // Task 1 assertion: the new fields exist, defaults (sentinels) reproduce the
  // known ec_37.a1 results computed via the plain Lfunc_init path.
  Lerror_t ec;
  Lfunc_t L = build_37(0.0, 0, "build/wt_cache_default", &ec);
  assert(!fatal_error(ec));
  assert(Lfunc_rank(L) == 1);
  // Cross-check the default-window zeros against LMFDB (bead 31c.8 validation):
  // EllipticCurve/Q/37/a positive_zeros lists 13 zeros (to t~19.81) at ~16 sig figs;
  // allow 2^-44 (~5.7e-14) for LMFDB's last-digit uncertainty. The enlarge run (task3)
  // reproduces these and extends past LMFDB's range, validated by overlap + |eps|=1.
  static const char *lmfdb_zeros[13] = {
    "5.003170014006659","6.870391216954432","8.014330807872879","9.933098353605352",
    "10.77513816254080","11.75732472284978","12.95838641388285","15.60385787320432",
    "16.19201741687448","17.14169364801487","18.06365420291071","18.78719562466392",
    "19.81482224536338"};
  for (int i = 0; i < 13; ++i) {
    arb_t ref; arb_init(ref);
    arb_set_str(ref, lmfdb_zeros[i], 300);
    arb_add_error_2exp_si(ref, -44);
    assert(arb_overlaps(Lfunc_zeros(L,0) + i, ref));
    arb_clear(ref);
  }
  printf("task1 ok (13 zeros vs LMFDB)\n");

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

  // ---- conv_support must be conductor-aware (Reviewer A blocker) ----
  // For conductor > ~1e4 (M0 = ceil(sqrt(N)/100) > 1) the small-coefficient path
  // (finish_convolves) fills the convolution down to bucket calc_m(1) =
  // round(log(1/sqrt(N))/(2pi/B)); once sqrt(N) > 100 that is BELOW the old
  // floor(log(0.01)/(2pi/B)) bound, so a conductor-independent conv_support sized
  // fft_N too small, the cyclic convolution aliased, and the certified zeros
  // silently stopped containing the truth -- with no error raised. Conductor 37
  // (M0=1, used by task1-3) never exercised this. Here the default window is
  // correct (it never aliases at any conductor) and the enlarged window's trusted
  // zeros must overlap it; on the buggy conv_support the enlarge aliases and the
  // overlap fails. fft_N goes 2048 (default) -> 4096 (enlarged at this conductor).
  {
    Lerror_t ecd;
    Lfunc_t Ld = build_dir(0.0, 0, "build/wt_cache_dir_def", &ecd);          // reference
    assert(!fatal_error(ecd));
    { arb_t m; arb_init(m); acb_abs(m, Lfunc_sign(Ld), 100);
      arb_sub_ui(m, m, 1, 100); assert(arb_contains_zero(m)); arb_clear(m); }
    int nd = 0; while (nd < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(Ld,0)+nd)) nd++;
    assert(nd > 0);

    Lerror_t ece;
    Lfunc_t Lde = build_dir(48.0, (uint64_t)1<<18, "build/wt_cache_dir_big", &ece);  // enlarge
    assert(!fatal_error(ece));
    { arb_t m; arb_init(m); acb_abs(m, Lfunc_sign(Lde), 100);
      arb_sub_ui(m, m, 1, 100); assert(arb_contains_zero(m)); arb_clear(m); }
    int ne = 0; while (ne < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(Lde,0)+ne)) ne++;
    assert(ne > nd);                                  // the wider window finds strictly more
    for (int i = 0; i < nd; ++i)                      // and every default zero reappears
      assert(arb_overlaps(Lfunc_zeros(Lde,0)+i, Lfunc_zeros(Ld,0)+i));
    Lfunc_clear(Ld); Lfunc_clear(Lde);
  }
  printf("highcond ok\n");

  // ---- df_zero static (pi*A)^d cache must not go stale across windows (bead 1yx) ----
  // df_zero (src/rank.c) recomputed a_pi_d[i] = (pi*A)^i only when the working
  // precision INCREASED -- keyed on prec, not on A. A = fft_NN/B varies per window,
  // so computing an enlarged object (higher precision, larger A) before a default
  // object of the SAME degree in one process left the default reading the enlarged
  // object's (pi*A)^d: its Lfunc_Taylor came back scaled by A_enlarge/A_default
  // (0.408 vs 0.306 here), a certified ball NOT containing the truth, raised with
  // only a non-fatal warning. The other tests miss it because none computes a
  // lower-precision window after a higher-precision one. (The companion turing.c
  // stale-cache half of 1yx -- a spurious, non-fatal ERR_RH_ERROR on the second
  // object -- is owned by branch fix/turing-static-cache and is tolerated here.)
  {
    Lerror_t ece, ecd;
    Lfunc_t Le2 = build_37(48.0, (uint64_t)1<<18, "build/wt_cache_1yx_big", &ece); // enlarge FIRST
    assert(!fatal_error(ece));
    Lfunc_t Ld2 = build_37(0.0, 0, "build/wt_cache_1yx_def", &ecd);                 // default SECOND
    assert(!fatal_error(ecd));
    arb_t bsd; arb_init(bsd);
    arb_set_str(bsd, "0.305999773834052301820483683321676474452637774590771998534541832481", 300);
    arb_add_error_2exp_si(bsd, -100); // generous: the bug is a ~33% scale error, not ulps
    assert(arb_overlaps(Lfunc_Taylor(Le2), bsd));  // enlarged object's leading Taylor is correct
    assert(arb_overlaps(Lfunc_Taylor(Ld2), bsd));  // default-after-enlarge correct (fails on broken df_zero)
    arb_clear(bsd);
    Lfunc_clear(Le2); Lfunc_clear(Ld2);
  }
  printf("multiwindow-taylor ok\n");

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

  // bead 31c.5: a valid shrink (H=16) returns a correct prefix of the default zeros.
  // H=16 (B=128) clears every lower-bound guard; its zeros reach only ~1.5*16=24,
  // a strict subset of the default (H=32, reach ~48). The low zeros must agree.
  {
    Lerror_t ecs;
    Lfunc_t Ls = build_37(16.0, 0, "build/wt_cache_shrink", &ecs);
    assert(!fatal_error(ecs));
    assert(Lfunc_rank(Ls) == 1);
    // |eps| = 1
    { arb_t m; arb_init(m); acb_abs(m, Lfunc_sign(Ls), 100);
      arb_sub_ui(m, m, 1, 100); assert(arb_contains_zero(m)); arb_clear(m); }
    int nshr = 0; while (nshr < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(Ls,0)+nshr)) nshr++;
    int ndef = 0; while (ndef < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(L,0)+ndef)) ndef++;
    assert(nshr > 0 && nshr < ndef);                 // strictly fewer zeros than default
    for (int i = 0; i < nshr; ++i)                   // and a correct prefix (overlap)
      assert(arb_overlaps(Lfunc_zeros(Ls,0)+i, Lfunc_zeros(L,0)+i));
    // the largest shrunk zero is within the H=16 reach (well under the default's ~48)
    arb_t lim; arb_init(lim); arb_set_d(lim, 30.0);
    assert(arb_lt(Lfunc_zeros(Ls,0)+(nshr-1), lim));
    arb_clear(lim);
    Lfunc_clear(Ls);
  }
  printf("shrink ok\n");

  // ---- special_value + Taylor are window-invariant AND match known values ----
  // The existing enlarge test (task3) only checks ZEROS overlap between default and
  // enlarged. Lfunc_special_value and Lfunc_Taylor exercise different windowed code:
  // both upsample off u_values_off using L->A = fft_NN/B, which scales with the window
  // (default fft_NN=2^16 vs enlarged 2^17 here, with a smaller 1/B). The L-value at a
  // fixed point and the central Taylor coefficient are mathematical invariants of the
  // L-function, so they MUST agree across windows -- and must contain the independently
  // known constants (LMFDB / examples/ec_37.a1.cpp):
  //   L(1.5) = 0.18396547525832984973...  (algebraic arg 1.5; analytic 1.0)
  //   L'(1/2)/1! = 0.305999773834052301820483683321676474452637774590772 (BSD).
  // If the window-dependent A (or anything it feeds) were wrong at the enlarged window,
  // the value would shift off the golden constant and/or stop overlapping the default.
  {
    Lerror_t ecd, ece;
    Lfunc_t Lsv_def = build_37(0.0, 0, "build/wt_cache_sv_def", &ecd);
    assert(!fatal_error(ecd));
    Lfunc_t Lsv_enl = build_37(48.0, (uint64_t)1<<18, "build/wt_cache_sv_enl", &ece);
    assert(!fatal_error(ece));

    acb_t vdef, venl; acb_init(vdef); acb_init(venl);
    assert(!fatal_error(Lfunc_special_value(vdef, Lsv_def, 1.5, 0.0)));
    assert(!fatal_error(Lfunc_special_value(venl, Lsv_enl, 1.5, 0.0)));

    // Each windowed value contains the known L(1.5).
    acb_t lref; acb_init(lref);
    arb_set_str(acb_realref(lref), "0.18396547525832984973", 300);
    arb_zero(acb_imagref(lref));
    arb_add_error_2exp_si(acb_realref(lref), -50);
    arb_add_error_2exp_si(acb_imagref(lref), -50);
    assert(acb_overlaps(vdef, lref));
    assert(acb_overlaps(venl, lref));
    // ... and they contain each other (window invariance of the L-value).
    assert(acb_overlaps(vdef, venl));
    acb_clear(lref);

    // The two computes are genuinely distinct objects, not the same handle: the
    // enlarged window's A differs, so its certified radius differs from the
    // default's (non-vacuous overlap). (default re-radius ~2.9e-39, enlarged ~8e-41.)
    double rd = mag_get_d(arb_radref(acb_realref(vdef)));
    double re = mag_get_d(arb_radref(acb_realref(venl)));
    assert(rd > 0.0 && re > 0.0 && rd != re);
    acb_clear(vdef); acb_clear(venl);

    // Taylor: central-point invariant, contains BSD, and overlaps across windows.
    arb_t bsd; arb_init(bsd);
    arb_set_str(bsd, "0.305999773834052301820483683321676474452637774590771998534541832481", 300);
    arb_add_error_2exp_si(bsd, -100);
    assert(arb_overlaps(Lfunc_Taylor(Lsv_def), bsd));
    assert(arb_overlaps(Lfunc_Taylor(Lsv_enl), bsd));
    assert(arb_overlaps(Lfunc_Taylor(Lsv_def), Lfunc_Taylor(Lsv_enl)));
    arb_clear(bsd);

    Lfunc_clear(Lsv_def); Lfunc_clear(Lsv_enl);
  }
  printf("special-value-window ok\n");

  // ---- equivalent (mus, normalisation) at a NON-default window => equal output ----
  // CLAUDE.md: for a weight-2 EC, {mus=[0,1], norm=0.5} and {mus=[0.5,1.5], norm=0}
  // describe the SAME analytic L-function (provided the Euler factors are given in
  // the matching normalisation -- see cb_anal). The window geometry must depend only
  // on the analytic data (L->mus = mus + norm = {0.5,1.5} both ways, and 1/B from H),
  // never on how the caller split it. So at H=48 the rank, every zero, and the central
  // Taylor coefficient must coincide. Non-vacuous: the two inputs drive different code
  // (coeff.c normalises by p^{-m*norm}; the mus split differs) -- supplying the SAME
  // algebraic coeffs with norm=0 instead yields a DIFFERENT L-function (rank 0 here),
  // so a regression that mixed the split into the geometry would break the overlap.
  {
    static double mus_a[2] = {0.0, 1.0};       // norm 0.5  -> analytic mus {0.5,1.5}
    static double mus_b[2] = {0.5, 1.5};       // norm 0.0  -> analytic mus {0.5,1.5}
    Lerror_t eca, ecb;
    Lfunc_t La = build_37_full(48.0, (uint64_t)1<<18, 0.5, mus_a, DK, DK, cb,
                               "build/wt_cache_eqv_a", &eca);
    assert(!fatal_error(eca));
    Lfunc_t Lb = build_37_full(48.0, (uint64_t)1<<18, 0.0, mus_b, DK, DK, cb_anal,
                               "build/wt_cache_eqv_b", &ecb);
    assert(!fatal_error(ecb));

    assert(Lfunc_rank(La) == 1);
    assert(Lfunc_rank(Lb) == 1);
    int na = wt_count_zeros(La, 0), nb = wt_count_zeros(Lb, 0);
    assert(na > 0 && na == nb);                 // same number of zeros found
    for (int i = 0; i < na; ++i)
      assert(arb_overlaps(Lfunc_zeros(La,0)+i, Lfunc_zeros(Lb,0)+i));
    // central Taylor coefficient (analytic central point is 1/2 in both splits)
    assert(arb_overlaps(Lfunc_Taylor(La), Lfunc_Taylor(Lb)));
    // |eps| = 1 for both
    for (Lfunc_t Lx : {La, Lb}) {
      arb_t m; arb_init(m); acb_abs(m, Lfunc_sign(Lx), 100);
      arb_sub_ui(m, m, 1, 100); assert(arb_contains_zero(m)); arb_clear(m);
    }
    Lfunc_clear(La); Lfunc_clear(Lb);
  }
  printf("equiv-normalisation ok\n");

  // ---- self_dual declared YES vs computed (DK), and side 0 vs side 1, at H=48 ----
  // ec_37.a1 IS self-dual. Declaring self_dual=YES makes compute.c skip find_zeros
  // for the dual (side 1); leaving it DK computes both sides. The certified side-0
  // output (rank, zeros, Taylor) must be identical whether or not the dual is also
  // computed -- the declaration only saves work, it must not change the answer. And
  // because the object is self-dual, side 1 (built from res[(-nn)%fft_NN] in copy())
  // must reproduce side 0 exactly. Both are strong cross-checks at a non-default
  // window. Non-vacuous: side 1 is computed from a different index path than side 0,
  // and the YES build takes a different control-flow branch.
  {
    static double mus[2] = {0,1};
    Lerror_t edk, eyes;
    Lfunc_t Ldk = build_37_full(48.0, (uint64_t)1<<18, 0.5, mus, DK, DK, cb,
                                "build/wt_cache_sd_dk", &edk);
    assert(!fatal_error(edk));
    Lfunc_t Lyes = build_37_full(48.0, (uint64_t)1<<18, 0.5, mus, YES, DK, cb,
                                 "build/wt_cache_sd_yes", &eyes);
    assert(!fatal_error(eyes));

    // YES vs DK agree on side 0.
    assert(Lfunc_rank(Lyes) == Lfunc_rank(Ldk));
    int ndk = wt_count_zeros(Ldk, 0), nyes = wt_count_zeros(Lyes, 0);
    assert(ndk > 0 && ndk == nyes);
    for (int i = 0; i < ndk; ++i)
      assert(arb_overlaps(Lfunc_zeros(Lyes,0)+i, Lfunc_zeros(Ldk,0)+i));
    assert(arb_overlaps(Lfunc_Taylor(Lyes), Lfunc_Taylor(Ldk)));

    // self_dual=YES skipped the dual: side 1 is empty. DK computed it: side 1 is not.
    assert(wt_count_zeros(Lyes, 1) == 0);
    int ndk1 = wt_count_zeros(Ldk, 1);
    assert(ndk1 == ndk);
    // side 0 and side 1 of the self-dual object coincide zero-for-zero.
    for (int i = 0; i < ndk; ++i)
      assert(arb_overlaps(Lfunc_zeros(Ldk,0)+i, Lfunc_zeros(Ldk,1)+i));

    Lfunc_clear(Ldk); Lfunc_clear(Lyes);
  }
  printf("self-dual-window ok\n");

  // ---- supplied rank vs computed rank at H=48 => same zeros / Taylor ----
  // Telling the library rank=1 up front (its true rank) must not change the certified
  // output versus letting it determine the rank (DK). The supplied value only short-
  // circuits do_rank; the zeros and the leading Taylor coefficient must be unchanged.
  // Non-vacuous: rank feeds the central-derivative bookkeeping (compute.c divides
  // L_d by i! up to rank), so a regression coupling rank into the wrong place would
  // shift Lfunc_Taylor off the DK value and off BSD.
  {
    static double mus[2] = {0,1};
    Lerror_t edk, er1;
    Lfunc_t Ldk = build_37_full(48.0, (uint64_t)1<<18, 0.5, mus, DK, DK, cb,
                                "build/wt_cache_rk_dk", &edk);
    assert(!fatal_error(edk));
    Lfunc_t Lr1 = build_37_full(48.0, (uint64_t)1<<18, 0.5, mus, DK, 1, cb,
                                "build/wt_cache_rk_1", &er1);
    assert(!fatal_error(er1));

    assert(Lfunc_rank(Lr1) == 1 && Lfunc_rank(Ldk) == 1);
    int ndk = wt_count_zeros(Ldk, 0), nr1 = wt_count_zeros(Lr1, 0);
    assert(ndk > 0 && ndk == nr1);
    for (int i = 0; i < ndk; ++i)
      assert(arb_overlaps(Lfunc_zeros(Lr1,0)+i, Lfunc_zeros(Ldk,0)+i));
    assert(arb_overlaps(Lfunc_Taylor(Lr1), Lfunc_Taylor(Ldk)));
    arb_t bsd; arb_init(bsd);
    arb_set_str(bsd, "0.305999773834052301820483683321676474452637774590771998534541832481", 300);
    arb_add_error_2exp_si(bsd, -100);
    assert(arb_overlaps(Lfunc_Taylor(Lr1), bsd));
    arb_clear(bsd);
    Lfunc_clear(Ldk); Lfunc_clear(Lr1);
  }
  printf("supplied-rank-window ok\n");

  Lfunc_clear(L);
  return 0;
}
