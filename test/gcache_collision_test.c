/*
  Regression test for the G-cache filename collision bug.

  compute_g() caches the gamma-factor (G) data to a file named after the
  L-function's mu's.  The filename must encode *all* the mu's, so that two
  L-functions with different mu-multisets (in particular different degrees)
  never share a cache file.

  The bug: the filename was built with a self-overlapping
      sprintf(fname1, "%s_%.1f", fname1, mu)
  whose undefined behaviour collapsed the name down to only the last mu.
  As a result a degree-2 L-function with analytic mu's [0.5, 1.5] and a
  degree-4 L-function with analytic mu's [0.5, 0.5, 1.5, 1.5] both mapped to
  the cache file "g_1.5".  Whichever ran first wrote the file; the second
  silently read the wrong G data and reported a wrong nmax (hence a wrong
  number of Euler factors / wrong bound).

  This test computes nmax for the degree-2 function twice:
    - once in a directory already populated by the degree-4 function
      (the collision scenario), and
    - once in a pristine directory (the correct reference value).
  A correct implementation gives the same nmax both times.
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

// Compute Lfunc_nmax for the given gamma data, using cache_dir for the
// G cache.  Uses conductor 1 so nmax reflects the (conductor-independent)
// exp(2 pi (hi_i + 1/2) / B) term directly.
static uint64_t nmax_in(uint64_t degree, double normalisation,
                        const double *mus, const char *cache_dir) {
  Lparams_t Lp = {0}; // zero-init: future Lparams_t fields default safely
  Lp.degree = degree;
  Lp.conductor = 1;
  Lp.normalisation = normalisation;
  Lp.mus = (double *)mus;          // init_advanced only reads from this
  Lp.target_prec = DEFAULT_TARGET_PREC;  // required for the cache path
  Lp.rank = DK;
  Lp.self_dual = DK;
  Lp.cache_dir = (char *)cache_dir;
  Lp.gprec = 0;                    // required for the cache path
  Lp.wprec = 0;
  Lp.max_t = 0;                    // default output window (64/degree)
  Lp.max_fft_NN = 0;               // default transform-size cap (1<<16)

  Lerror_t ecode = ERR_SUCCESS;
  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  assert(!fatal_error(ecode));
  uint64_t M = Lfunc_nmax(L);
  Lfunc_clear(L);
  return M;
}

int main(void) {
  // analytic mu's are mus[i] + normalisation:
  //   degree 4: [0,0,1,1] + 0.5 = [0.5, 0.5, 1.5, 1.5]
  //   degree 2: [0,1]     + 0.5 = [0.5, 1.5]
  // -> both share the largest mu 1.5, so the buggy name collapses both to g_1.5.
  double mus4[4] = {0.0, 0.0, 1.0, 1.0};
  double mus2[2] = {0.0, 1.0};

  char shared_dir[] = "./gcache_collide_XXXXXX";
  char ref_dir[]    = "./gcache_ref_XXXXXX";
  assert(mkdtemp(shared_dir) != NULL);
  assert(mkdtemp(ref_dir)    != NULL);

  // Populate the shared dir with the degree-4 cache first.
  uint64_t m4 = nmax_in(4, 0.5, mus4, shared_dir);
  // Now the degree-2 function in the same dir: must NOT pick up g-data
  // that belongs to the degree-4 function.
  uint64_t m2_shared = nmax_in(2, 0.5, mus2, shared_dir);
  // The trusted degree-2 value, computed with no pre-existing cache.
  uint64_t m2_ref = nmax_in(2, 0.5, mus2, ref_dir);

  printf("nmax degree-4 [0.5,0.5,1.5,1.5] = %" PRIu64 "\n", m4);
  printf("nmax degree-2 [0.5,1.5] (shared cache dir) = %" PRIu64 "\n", m2_shared);
  printf("nmax degree-2 [0.5,1.5] (pristine cache dir) = %" PRIu64 "\n", m2_ref);

  // Clean up before asserting, so the temp dirs are removed on both the pass
  // and the (abort-on-assert) failure path.  The printed nmax values above are
  // enough to diagnose a failure.
  cleanup_dir(shared_dir);
  cleanup_dir(ref_dir);

  // Sanity: the two L-functions genuinely differ, so the test is meaningful.
  assert(m4 != m2_ref);
  // Regression: the degree-2 nmax must be independent of whether a degree-4
  // cache file is sitting in the directory.
  assert(m2_shared == m2_ref);

  printf("PASS: G cache does not collide across distinct L-functions\n");
  return 0;
}
