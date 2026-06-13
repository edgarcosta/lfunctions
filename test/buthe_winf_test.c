/*
   Acceptance test for buthe_winf_integral (src/buthe_winf.c).

   Checks the rigorous Arb enclosure of the per-Gamma-factor archimedean
   integral

     I(b,h,mu) = 2 * \int_0^\infty
       e^{-(1/2+mu) t} / (1 - e^{-2 t})
       * ( b/pi - sin(b t)/(pi t cosh(h t / 2)) )  dt

   against 22-significant-figure references (true I(b,h,mu) values).

   For each case we assert that the returned ball CONTAINS the reference
   and that its radius is < 1e-10.  Failure is signalled via assert.
*/

#include <assert.h>
#include <flint/arb.h>
#include "glfunc_internals.h"

struct testcase {
  double h;
  double b;
  double mu;
  const char *ref;
};

int main(void)
{
  const slong prec = 200;

  // 22 significant figures (true I(b,h,mu) values).
  struct testcase cases[] = {
    {8.0, 32.0,       0.0, "60.94706827230718282252"},
    {8.0, 32.0,       1.0, "29.44706600870153811482"},
    {8.0, 16.0,       0.0, "26.94969404613122873693"},
    {8.0, 16.0,       1.0, "11.44848190585916734350"},
    {8.0, 64.0 / 9.0, 0.0, "10.35318974619509695387"},

    {6.0, 32.0,       0.0, "60.91145311730560447573"},
    {6.0, 16.0,       1.0, "11.37141615129653364139"},
    {6.0, 64.0 / 9.0, 0.5, "5.290851192028140134374"},

    {5.0, 32.0,       0.5, "38.08850624782670443630"},
    {5.0, 16.0,       0.0, "26.84231945214516413477"},
    {5.0, 64.0 / 9.0, 1.0, "3.470991475273731660225"},

    {4.0, 32.0,       1.0, "29.38630512947141821112"},
    {4.0, 16.0,       0.5, "15.54135979908765070632"},
    {4.0, 64.0 / 9.0, 0.0, "10.02759007442487484771"},

    {3.0, 32.0,       0.0, "60.87755526878309439761"},
    {3.0, 16.0,       1.0, "11.30118455338333089650"},
    {3.0, 64.0 / 9.0, 0.5, "5.113020782373083911300"},
  };
  const size_t n = sizeof(cases) / sizeof(cases[0]);

  arb_t b, h, got, ref, rad, tol;
  arb_init(b);
  arb_init(h);
  arb_init(got);
  arb_init(ref);
  arb_init(rad);
  arb_init(tol);

  // tol = 1e-10
  arb_set_str(tol, "1e-10", prec);

  int ok = 1;
  for (size_t i = 0; i < n; i++) {
    arb_set_d(b, cases[i].b);
    arb_set_d(h, cases[i].h);
    buthe_winf_integral(got, b, h, cases[i].mu, prec);

    if (arb_set_str(ref, cases[i].ref, prec) != 0) {
      flint_printf("FAIL: could not parse reference '%s'\n", cases[i].ref);
      ok = 0;
      continue;
    }

    arb_get_rad_arb(rad, got);

    int contains = arb_contains(got, ref);
    int tight = arb_lt(rad, tol);

    flint_printf("h=%.1f b=%.10g mu=%.1f : ",
                 cases[i].h, cases[i].b, cases[i].mu);
    flint_printf("got = "); arb_printd(got, 20);
    flint_printf("  rad = "); arb_printd(rad, 6);
    flint_printf("  contains=%d tight=%d\n", contains, tight);

    if (!contains) {
      flint_printf("  FAIL: enclosure does not contain reference %s\n",
                   cases[i].ref);
      ok = 0;
    }
    if (!tight) {
      flint_printf("  FAIL: radius not < 1e-10\n");
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
