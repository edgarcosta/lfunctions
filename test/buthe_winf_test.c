/*
   Acceptance test for buthe_winf_integral (src/buthe_winf.c).

   Checks the rigorous Arb enclosure of the per-Gamma-factor archimedean
   integral

     I(b,h,mu) = 2 * \int_0^\infty
       e^{-(1/2+mu) t} / (1 - e^{-2 t})
       * ( b/pi - sin(b t)/(pi t cosh(h t / 2)) )  dt

   against 22-significant-figure references (true I(b,h=8,mu) values) for h = 8.

   For each case we assert that the returned ball CONTAINS the reference
   and that its radius is < 1e-9.  Failure is signalled via assert.
*/

#include <assert.h>
#include <flint/arb.h>
#include "glfunc_internals.h"

struct testcase {
  double b;
  double mu;
  const char *ref;
};

int main(void)
{
  const slong prec = 200;

  // h = 8 throughout (matches gp/buthe_ints.gp).
  const double h_val = 8.0;

  // 22 significant figures (true I(b,h=8,mu) values).
  struct testcase cases[] = {
    {32.0,       0.0, "60.94706827230718282252"},
    {32.0,       1.0, "29.44706600870153811482"},
    {16.0,       0.0, "26.94969404613122873693"},
    {16.0,       1.0, "11.44848190585916734350"},
    {64.0 / 9.0, 0.0, "10.35318974619509695387"},
  };
  const size_t n = sizeof(cases) / sizeof(cases[0]);

  arb_t b, h, got, ref, rad, tol;
  arb_init(b);
  arb_init(h);
  arb_init(got);
  arb_init(ref);
  arb_init(rad);
  arb_init(tol);

  // tol = 1e-9
  arb_set_str(tol, "1e-9", prec);

  arb_set_d(h, h_val);

  int ok = 1;
  for (size_t i = 0; i < n; i++) {
    arb_set_d(b, cases[i].b);
    buthe_winf_integral(got, b, h, cases[i].mu, prec);

    if (arb_set_str(ref, cases[i].ref, prec) != 0) {
      flint_printf("FAIL: could not parse reference '%s'\n", cases[i].ref);
      ok = 0;
      continue;
    }

    arb_get_rad_arb(rad, got);

    int contains = arb_contains(got, ref);
    int tight = arb_lt(rad, tol);

    flint_printf("b=%.10g mu=%.1f : ", cases[i].b, cases[i].mu);
    flint_printf("got = "); arb_printd(got, 20);
    flint_printf("  rad = "); arb_printd(rad, 6);
    flint_printf("  contains=%d tight=%d\n", contains, tight);

    if (!contains) {
      flint_printf("  FAIL: enclosure does not contain reference %s\n",
                   cases[i].ref);
      ok = 0;
    }
    if (!tight) {
      flint_printf("  FAIL: radius not < 1e-9\n");
      ok = 0;
    }

    assert(contains);
    assert(tight);
  }

  arb_clear(b);
  arb_clear(h);
  arb_clear(got);
  arb_clear(ref);
  arb_clear(rad);
  arb_clear(tol);

  if (!ok) {
    flint_printf("buthe_winf_test: FAILED\n");
    return 1;
  }
  flint_printf("buthe_winf_test: all %zu cases passed\n", n);
  return 0;
}
