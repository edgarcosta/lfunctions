/*
  Regression test for the G-cache staleness footgun (bead lfunctions-fe3).

  compute_g() caches the gamma-factor (G) data to a file named after the
  L-function's mu's.  The cached grid extent (low_i, hi_i, max_K) and the
  precision it was computed at (gprec) are a function of (degree, mus, gprec);
  hi_i and max_K grow monotonically with gprec.  The cache is consulted
  whenever the *input* gprec is 0 (default mode), but the *effective* gprec is
  then derived as max(target+gamma+EXTRA_BITS, wprec).  A caller who raises
  wprec (as one does for a large conductor) therefore needs a higher-precision
  grid than a default run wrote.

  The bug: read_gfile() validated nothing.  A cheap cache file written by a
  default (wprec=0) run was silently reused by a later run with an elevated
  wprec, even though its gprec -- hence hi_i / max_K -- was too small.  The
  result was a silently-truncated grid and wrong numbers, with no error.  The
  documented workaround was "rm -f g_* before running".

  The fix: the cache file carries a self-describing header (magic, version,
  degree, the sorted mus as exact halved integers, and the gprec it was
  computed at).  On load the library verifies the file is USABLE for the
  current request -- same degree and mus, and cached gprec >= required gprec --
  and, when it is not, recomputes and overwrites the stale file instead of
  returning wrong numbers.  (A file whose header is valid but whose body is
  truncated/corrupt is a genuinely broken file and stays a fatal ERR_G_INFILE.)

  We probe at init time via Lfunc_nmax, which reflects hi_i and so reveals an
  undersized grid, exactly as test/gcache_collision_test.c does.
*/

#include <assert.h>
#include <dirent.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "glfunc.h"

// Remove a flat directory (the G cache writes only regular files into it).
static void cleanup_dir(const char *dir) {
  DIR *d = opendir(dir);
  if (d) {
    struct dirent *e;
    while ((e = readdir(d)) != NULL) {
      if (!strcmp(e->d_name, ".") || !strcmp(e->d_name, "..")) continue;
      char path[1200];
      snprintf(path, sizeof(path), "%s/%s", dir, e->d_name);
      remove(path);
    }
    closedir(d);
  }
  rmdir(dir);
}

// Initialise an L-function in default cache mode (gprec=0) with the given
// working precision and cache directory.  Returns whether init reported a
// fatal error; when it did not, *out_nmax receives Lfunc_nmax (which reflects
// the G-grid extent).  conductor 1 keeps nmax conductor-independent.
static bool init_is_fatal(uint64_t degree, double normalisation,
                          const double *mus, int64_t wprec,
                          const char *cache_dir, uint64_t *out_nmax) {
  Lparams_t Lp;
  Lp.degree = degree;
  Lp.conductor = 1;
  Lp.normalisation = normalisation;
  Lp.mus = (double *)mus;                // init_advanced only reads from this
  Lp.target_prec = DEFAULT_TARGET_PREC;  // required for the cache path
  Lp.rank = DK;
  Lp.self_dual = DK;
  Lp.cache_dir = (char *)cache_dir;
  Lp.gprec = 0;                          // required for the cache path
  Lp.wprec = wprec;

  Lerror_t ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  bool fatal = fatal_error(ecode);
  if (out_nmax) *out_nmax = (!fatal && L) ? Lfunc_nmax(L) : 0;
  if (L) Lfunc_clear(L);               // init returns NULL on fatal error
  return fatal;
}

int main(void) {
  // degree 2, analytic mu's [0.5, 1.5] (= [0,1] + 0.5), as for a weight-2
  // elliptic curve.
  double mus2[2] = {0.0, 1.0};
  const double norm = 0.5;
  // A working precision well above the default (~target+EXTRA_BITS bits) so
  // the effective gprec -- and hence the required grid -- is strictly larger
  // than a default run writes.
  const int64_t HI_WPREC = 512;

  char dir_def[]   = "./gcache_def_XXXXXX";    // pristine, default precision
  char dir_hi[]    = "./gcache_hi_XXXXXX";     // pristine, elevated wprec
  char dir_stale[] = "./gcache_stale_XXXXXX";  // default then elevated (reuse)
  assert(mkdtemp(dir_def)   != NULL);
  assert(mkdtemp(dir_hi)    != NULL);
  assert(mkdtemp(dir_stale) != NULL);

  // Reference nmax for a default run and for an elevated-wprec run, each in a
  // pristine cache directory.
  uint64_t nmax_default = 0, nmax_hi = 0;
  bool fatal_def = init_is_fatal(2, norm, mus2, 0,        dir_def, &nmax_default);
  bool fatal_hi  = init_is_fatal(2, norm, mus2, HI_WPREC, dir_hi,  &nmax_hi);

  // Populate dir_stale with the cheap default cache, then issue the
  // elevated-wprec request in the SAME directory: it must NOT silently reuse
  // the undersized cache.
  uint64_t nmax_seed = 0, nmax_stale = 0;
  bool fatal_seed  = init_is_fatal(2, norm, mus2, 0,        dir_stale, &nmax_seed);
  bool fatal_stale = init_is_fatal(2, norm, mus2, HI_WPREC, dir_stale, &nmax_stale);

  printf("nmax default (wprec=0)         = %" PRIu64 " (fatal=%d)\n", nmax_default, fatal_def);
  printf("nmax elevated (wprec=%" PRId64 ")     = %" PRIu64 " (fatal=%d)\n", HI_WPREC, nmax_hi, fatal_hi);
  printf("nmax stale-reuse (elevated)    = %" PRIu64 " (fatal=%d)\n", nmax_stale, fatal_stale);

  cleanup_dir(dir_def);
  cleanup_dir(dir_hi);
  cleanup_dir(dir_stale);

  // The pristine runs must both succeed.
  assert(!fatal_def);
  assert(!fatal_hi);
  assert(!fatal_seed);

  // Sanity: the elevated-wprec run genuinely needs a larger grid than the
  // default run, so reusing the default cache really would be wrong.  If this
  // ever fails the test has stopped exercising the bug.
  assert(nmax_hi > nmax_default);

  // Regression: the elevated request issued against a directory holding only
  // the cheap default cache must NOT silently reuse it.  The insufficient file
  // is recomputed and overwritten, yielding the same deep grid (nmax) as a
  // pristine elevated run -- not a fatal error, and not the stale shallow grid.
  assert(!fatal_stale);
  assert(nmax_stale == nmax_hi);

  // Sufficiency, not equality: a grid written at HIGH precision must satisfy a
  // later DEFAULT request (cached gprec >= required gprec), so it is reused
  // rather than rejected.  Guards the fix against over-rejecting usable caches.
  char dir_over[] = "./gcache_over_XXXXXX";
  assert(mkdtemp(dir_over) != NULL);
  uint64_t nmax_seed_hi = 0, nmax_reuse_lo = 0;
  bool fatal_seed_hi  = init_is_fatal(2, norm, mus2, HI_WPREC, dir_over, &nmax_seed_hi);
  bool fatal_reuse_lo = init_is_fatal(2, norm, mus2, 0,        dir_over, &nmax_reuse_lo);
  cleanup_dir(dir_over);
  assert(!fatal_seed_hi);
  assert(!fatal_reuse_lo);      // over-sufficient cache reused, not rejected

  // No happy-path regression: a matching cache is reused and gives the same
  // answer as a pristine run.
  char dir_same[] = "./gcache_same_XXXXXX";
  assert(mkdtemp(dir_same) != NULL);
  uint64_t nmax_w1 = 0, nmax_w2 = 0;
  bool fatal_w1 = init_is_fatal(2, norm, mus2, 0, dir_same, &nmax_w1);
  bool fatal_w2 = init_is_fatal(2, norm, mus2, 0, dir_same, &nmax_w2);
  cleanup_dir(dir_same);
  assert(!fatal_w1);
  assert(!fatal_w2);
  assert(nmax_w1 == nmax_default);   // first (pristine) run matches reference
  assert(nmax_w2 == nmax_default);   // cache reuse gives the identical answer

  printf("PASS: stale G cache rejected; sufficient/matching caches reused\n");
  return 0;
}
