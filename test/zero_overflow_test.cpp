// Opt-in regression test: hitting the per-side zero cap MAX_ZEROS must fail
// LOUDLY and FATALLY with ERR_ZERO_OVERFLOW, never silently truncate.
//
// Before this fix, find_zeros (src/zeros.c) `return ecode` when the zeros array
// filled, indistinguishable from the legitimate end-of-range return. A truncated
// zero set then fed the RH check, silently invalidating it. The fix flags
// ERR_ZERO_OVERFLOW (fatal) at the three cap sites so a caller can never trust a
// truncated result.
//
// This drives a REAL degree-2 computation whose analyzed range packs MORE than
// MAX_ZEROS (2048) zeros, and asserts the result is fatal with ERR_ZERO_OVERFLOW
// set. It is HEAVY (fft_NN = 2^20, ~1.4M Euler factors, ~3 GB peak resident;
// the zero-finding alone isolates ~2048 zeros, so it runs for several minutes,
// longer on a busy machine), so it is DENY-LISTED from the default `make check`
// via CHECK_EXCLUDE in Makefile.in. Build and run it directly:
//
//     ./configure --enable-assert && make test
//     ./build/test/zero_overflow_test.exe
//
// Exit 0 (no assertion fired) = pass.
//
// The trigger object is the product of the quadratic characters mod 10007 and
// mod 10009 (conductor 100160063 ~ 1e8), a self-dual degree-2 L-function, run at
// an enlarged output window H = max_t = 384 with target_prec = 300. A large
// conductor boosts the zero DENSITY (~ (1/2pi) log N per unit height) without
// growing fft_NN as fast as raising max_t on a small conductor would, and the
// raised precision lets the working precision survive to the far (1.5*H) end of
// the analyzed range where the 2048nd zero lies. At these parameters find_zeros
// fills the array (count == MAX_ZEROS) before reaching the end of the Turing
// zone, so the cap, not legitimate completion, is what stops it.

#include <flint/acb_poly.h>
#include <flint/ulong_extras.h>
#include "glfunc.h"
#include <cassert>
#include <cstdio>
#include <cstdint>
#include <sys/stat.h>

// Good-prime Euler factor of chi_P * chi_Q at p: (1 - chi_P(p) x)(1 - chi_Q(p) x).
// A ramified chi(p) = 0 drops its factor (n_jacobi returns 0 on a shared factor).
static long P = 10007, Q = 10009;
static void dir_cb(acb_poly_t poly, uint64_t p, int, int64_t, void *) {
  acb_poly_zero(poly);
  long a = n_jacobi((slong)(p % P), P);   // chi_P(p) in {-1,0,1}
  long b = n_jacobi((slong)(p % Q), Q);   // chi_Q(p)
  acb_poly_set_coeff_si(poly, 0, 1);
  acb_poly_set_coeff_si(poly, 1, -(a + b));
  if (a != 0 && b != 0) acb_poly_set_coeff_si(poly, 2, a * b);
}

// Count the (contiguous, zero-terminated) zeros stored on a side.
static int count_zeros(Lfunc_t L, uint64_t side) {
  int n = 0;
  while (n < (int)MAX_ZEROS && !arb_is_zero(Lfunc_zeros(L, side) + n)) n++;
  return n;
}

int main() {
  static double mus[2] = {0, 1};
  Lparams_t Lp = {};
  Lp.degree = 2;
  Lp.conductor = (uint64_t)P * (uint64_t)Q; // 100160063 ~ 1e8
  Lp.normalisation = 0.0;
  Lp.mus = mus;
  Lp.target_prec = 300;                     // survive precision decay to the 1.5*H Turing end
  Lp.wprec = 0;
  Lp.gprec = 0;
  Lp.self_dual = YES;                        // self-dual: skip the dual side (halve the work)
  Lp.rank = DK;
  Lp.cache_dir = (char *)"build/zero_overflow_cache";
  Lp.max_t = 384.0;                          // H: want_fft_NN = next_pow2(2048*384) = 2^20
  Lp.max_fft_NN = (uint64_t)1 << 21;         // cap comfortably above the required 2^20

  mkdir("build/zero_overflow_cache", 0777);

  Lerror_t ec = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ec);
  assert(!fatal_error(ec));                  // the window/geometry itself is valid

  ec |= Lfunc_use_all_lpolys(L, dir_cb, NULL);
  assert(!fatal_error(ec));                  // supplying Euler factors must not fail

  ec |= Lfunc_compute(L);

  // The discriminating assertions: on the FIXED code the cap is reached and the
  // result is FATAL with ERR_ZERO_OVERFLOW set; on the BROKEN code (silent
  // `return ecode` at the cap sites) the same compute fills the array but returns
  // a NON-fatal code, so both of these fail.
  assert(fatal_error(ec));
  assert(ec & ERR_ZERO_OVERFLOW);

  // Corroboration: the array did fill to exactly the cap (this alone does NOT
  // distinguish fixed from broken -- the broken code also stores MAX_ZEROS zeros
  // -- so it is a sanity check, not the discriminator).
  int nz = count_zeros(L, 0);
  assert(nz == (int)MAX_ZEROS);

  Lfunc_clear(L);
  printf("zero_overflow ok (found %d zeros on side 0, ERR_ZERO_OVERFLOW fatal)\n", nz);
  return 0;
}
