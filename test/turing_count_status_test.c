// Unit test for the Turing zero-count interval classifier.
//
// The old RH check used the analytic LOWER bound of N(t0) to declare
// "too many zeros" when lo > zeros_found - 1.  That misclassifies every wide
// but otherwise consistent count interval with lo above zeros_found - 1.  The
// corrected test only declares "too many" when the UPPER bound is below the
// rigorous number of isolated zeros.

#include <assert.h>
#include <flint/arb.h>
#include "glfunc_internals.h"

static void set_interval_d(arb_t x, double mid, double rad)
{
  arb_t err;
  arb_init(err);
  arb_set_d(x, mid);
  arb_set_d(err, rad);
  arb_add_error(x, err);
  arb_clear(err);
}

int main(void)
{
  arb_t tcount;
  arb_init(tcount);

  // Exact count: N(t0) = zeros_found.
  arb_set_ui(tcount, 5);
  assert(turing_count_status(tcount, 5, DEFAULT_TARGET_PREC) == TURING_COUNT_CONFIRMED);

  // Wide but certifying: [4.4, 5.6] contains no integer >= 5 except 5.
  // The old lower-bound check saw lo > 4 and falsely called this "too many".
  set_interval_d(tcount, 5.0, 0.6);
  assert(turing_count_status(tcount, 5, DEFAULT_TARGET_PREC) == TURING_COUNT_CONFIRMED);

  // Wide and incomplete: [5.4, 6.4] is consistent with a missed sixth zero.
  // This is "too few", not "too many", despite lo > zeros_found - 1.
  set_interval_d(tcount, 5.9, 0.5);
  assert(turing_count_status(tcount, 5, DEFAULT_TARGET_PREC) == TURING_COUNT_TOO_FEW);

  // Impossible over-count: the analytic upper bound is below the found zeros.
  set_interval_d(tcount, 3.9, 0.5);
  assert(turing_count_status(tcount, 5, DEFAULT_TARGET_PREC) == TURING_COUNT_TOO_MANY);

  arb_clear(tcount);
  return 0;
}
