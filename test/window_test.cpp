// Regression test for the configurable output window (bead lfunctions-31c).
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"
#include <cassert>
#include <cstdio>
#include <map>
#include <vector>
#include <cstdint>
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
  *ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  if (fatal_error(*ecode)) return L;
  *ecode |= Lfunc_use_all_lpolys(L, cb, NULL);
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

  Lfunc_clear(L);
  return 0;
}
