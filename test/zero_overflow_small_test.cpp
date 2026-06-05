// Fast default-suite regression for zero-cap exhaustion.
//
// The heavy zero_overflow_test.cpp proves that the real MAX_ZEROS=2048 cap is
// reachable and fatal on a large computation. This test uses the internal test
// hook to lower the active cap on an ordinary ec_37.a1 run, so make check covers
// the same control-flow invariant cheaply: filling the zero storage is fatal and
// sets ERR_ZERO_OVERFLOW.

#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"
#include <cassert>
#include <cstdint>
#include <map>
#include <vector>
#include <sys/stat.h>
using std::map; using std::vector;

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

static int count_zeros_to_cap(Lfunc_t L, uint64_t side, int cap) {
  int n = 0;
  while (n < cap && !arb_is_zero(Lfunc_zeros(L, side) + n)) n++;
  return n;
}

int main() {
  static double mus[2] = {0, 1};
  Lparams_t Lp;
  Lparams_init(&Lp);
  Lp.degree = 2;
  Lp.conductor = 37;
  Lp.normalisation = 0.5;
  Lp.mus = mus;
  Lp.self_dual = YES;
  Lp.cache_dir = (char*)"build/zero_overflow_small_cache";
  mkdir("build/zero_overflow_small_cache", 0777);

  Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  assert(!fatal_error(ec));
  Lfunc_set_zero_cap_for_tests(L, 2);

  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  assert(!fatal_error(ec));

  ec |= Lfunc_compute(L);
  assert(fatal_error(ec));
  assert(ec & ERR_ZERO_OVERFLOW);
  assert(count_zeros_to_cap(L, 0, 2) == 2);

  Lfunc_clear(L);
  return 0;
}
