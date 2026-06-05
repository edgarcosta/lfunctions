// Init-only regression tests for configurable output-window geometry.
#include "glfunc.h"
#include "glfunc_internals.h"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>
using std::vector;

// Synthetic degree-r object (mus all 0, normalisation 0.5 => analytic mus all
// 0.5) at a chosen window/cap. B / one_over_B / fft_NN / A depend only on degree
// and H, so no Euler factors or compute are needed. A non-zero gprec bypasses
// the on-disk G-cache and keeps this init-only sweep fast.
static Lfunc_t geo_init(uint64_t degree, double max_t, uint64_t max_fft_NN, Lerror_t *ecode) {
  vector<double> mus(degree, 0.0);
  Lparams_t Lp;
  Lparams_init(&Lp);
  Lp.degree = degree;
  Lp.conductor = 11;
  Lp.normalisation = 0.5;
  Lp.mus = mus.data();
  Lp.gprec = 30;
  Lp.cache_dir = NULL;
  Lp.max_t = max_t;
  Lp.max_fft_NN = max_fft_NN;
  *ecode = ERR_SUCCESS;
  return Lfunc_init_advanced(&Lp, ecode); // copies mus
}

int main() {
  // Power-of-two geometry sweep. For degree r and H = 2^k/r, the geometry is
  // forced to exact powers of two:
  //   B = 2^(k+3)/r, one_over_B = r/2^(k+3), fft_NN = 2^(k+10), A = 128*r.
  for (uint64_t r = 2; r <= (uint64_t)MAX_DEGREE; ++r) {
    for (int k = 1; k <= 9; ++k) {
      double H = ldexp(1.0, k) / (double)r;
      Lerror_t ecg;
      Lfunc_t Lg = geo_init(r, H, (uint64_t)1 << 20, &ecg);
      if (fatal_error(ecg)) {
        assert(ecg & ERR_WINDOW_TOO_SMALL);
        if (Lg) Lfunc_clear(Lg);
        continue;
      }
      Lfunc *G = (Lfunc*)Lg;
      assert(G->fft_NN == ((uint64_t)1 << (k + 10)));
      assert(G->one_over_B == (double)r / ldexp(1.0, k + 3));
      assert(G->A == 128.0 * (double)r);
      double reach_out = (double)G->fft_NN / ((double)OUTPUT_RATIO * G->A);
      assert(reach_out == H);
      double reach_tur = ((double)(G->fft_NN / OUTPUT_RATIO) + (double)(G->fft_NN / TURING_RATIO)) / G->A;
      assert(fabs(reach_tur - 1.5 * H) <= 1e-12 * (1.5 * H));
      Lfunc_clear(Lg);
    }
  }
  printf("geometry-sweep ok\n");

  // Explicit H = 64/r reproduces the sentinel default geometry.
  for (uint64_t r = 2; r <= (uint64_t)MAX_DEGREE; ++r) {
    Lerror_t es, ee;
    Lfunc_t Ls = geo_init(r, 0.0, 0, &es);
    Lfunc_t Le = geo_init(r, 64.0 / (double)r, (uint64_t)1 << 20, &ee);
    assert(!fatal_error(es) && !fatal_error(ee));
    Lfunc *S = (Lfunc*)Ls, *E = (Lfunc*)Le;
    assert(S->fft_NN == ((uint64_t)1 << 16) && E->fft_NN == S->fft_NN);
    assert(E->one_over_B == S->one_over_B);
    assert(E->A == S->A && S->A == 128.0 * (double)r);
    Lfunc_clear(Ls); Lfunc_clear(Le);
  }
  printf("explicit-default ok\n");

  // H = 128/r needs fft_NN = 2^17: a 2^16 cap rejects, a 2^17 cap admits.
  for (uint64_t r = 2; r <= (uint64_t)MAX_DEGREE; ++r) {
    double H = 128.0 / (double)r;
    Lerror_t e1, e2;
    Lfunc_t L1 = geo_init(r, H, (uint64_t)1 << 16, &e1);
    assert(e1 & ERR_WINDOW_TOO_LARGE);
    if (L1) Lfunc_clear(L1);
    Lfunc_t L2 = geo_init(r, H, (uint64_t)1 << 17, &e2);
    assert(!fatal_error(e2));
    assert(((Lfunc*)L2)->fft_NN == ((uint64_t)1 << 17));
    Lfunc_clear(L2);
  }
  printf("cap-threshold ok\n");

  // Just above 64/r the required length rounds up to 2^17; just below it stays
  // at 2^16. In both non-default cases A = fft_NN/(8H), not exactly 128*r.
  for (uint64_t r = 2; r <= (uint64_t)MAX_DEGREE; ++r) {
    double A128 = 128.0 * (double)r;
    Lerror_t eu, eb;
    Lfunc_t Lu = geo_init(r, 64.0 / (double)r + 1.0 / 1024.0, (uint64_t)1 << 20, &eu);
    assert(!fatal_error(eu));
    assert(((Lfunc*)Lu)->fft_NN == ((uint64_t)1 << 17) && ((Lfunc*)Lu)->A > A128);
    Lfunc_clear(Lu);
    Lfunc_t Lb = geo_init(r, 64.0 / (double)r - 1.0 / 1024.0, (uint64_t)1 << 20, &eb);
    assert(!fatal_error(eb));
    assert(((Lfunc*)Lb)->fft_NN == ((uint64_t)1 << 16) && ((Lfunc*)Lb)->A > A128);
    Lfunc_clear(Lb);
  }
  printf("round-boundary ok\n");

  return 0;
}
