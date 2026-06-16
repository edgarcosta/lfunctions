/*
   buthe_winf.c

   Rigorous on-the-fly Arb computation of the per-Gamma-factor archimedean
   integral used by Buthe's RH-verification inequality:

     I(b,h,mu) = 2 * \int_0^\infty
        e^{-(1/2+mu) t} / (1 - e^{-2 t})
        * ( b/pi - sin(b t)/(pi t cosh(h t / 2)) )  dt.

   This replaces the role of the static precomputed table in gp/buthe_ints.gp
   (captured as src/../gp/buthe_ints.out): instead of a fixed table indexed by
   degree and a half-integer mu, we evaluate the integral for arbitrary (b,h,mu)
   and return a verified enclosure (a ball that provably contains the true
   value).

   Regularity at the origin.
   -------------------------
   Write the integrand as  g(t) = E(t) * A(t) * B(t)  with
       E(t) = e^{-(1/2+mu) t},
       A(t) = 1 / (1 - e^{-2 t})            (simple pole at t = 0),
       B(t) = b/pi - sin(b t)/(pi t cosh(h t / 2)).
   B has a *double* zero at t = 0, so A*B (hence g) is regular there, with
   g(0) = 0.  A ball straddling 0, however, encloses the pole of A, so g cannot
   be evaluated naively on such a ball.  We therefore split

       [0, infty) = [0, t0]  U  [t0, T]  U  [T, infty),   0 < t0 < pi/8,

   and treat each piece by a method that never divides by an interval
   containing 0:

   * [0, t0]:  g is analytic on the disc |t| <= R0 for any R0 < pi/8 (the
     nearest singularity is the pole of 1/cosh(h t/2) at t = i pi/8 when h = 8;
     more generally the nearest singularity to 0 is at distance pi/h from the
     pole of cosh, or pi from the pole of A, whichever is smaller).  We expand g
     in a Taylor series at 0 using Arb power-series arithmetic.  The factors
     1/t inside A and B are handled by the exact "shift right" of a series whose
     constant term is provably 0 (1 - e^{-2t} and the inner bracket of B both
     vanish at 0), so no 0/0 is ever formed.  We integrate the truncated series
     exactly and add a closed-form Cauchy tail bound for the discarded terms.

   * [t0, T]:  g is regular on the whole interval (t0 > 0), so we hand it to
     Arb's rigorous integrator acb_calc_integrate, which returns a verified
     enclosure.  The integrand certifies analyticity (returns a finite ball
     when the input ball avoids every pole, and an indeterminate value
     otherwise, prompting the integrator to subdivide).

   * [T, infty):  bounded in closed form.  Since |sin| <= 1 and cosh >= 1,
       |g(t)| <= e^{-(1/2+mu) t} (b/pi + 1/(pi t)) / (1 - e^{-2 t}),
     and with alpha = 1/2 + mu > 0,
       \int_T^infty |g| <= (1/(1-e^{-2T})) * (e^{-alpha T}/alpha) * (b/pi + 1/(pi T)).
     T is chosen large enough that this bound is far below the target accuracy;
     correctness holds for any T > t0 because the bound is added as a symmetric
     error term.

   Finally the three contributions are summed and multiplied by 2.

   Accuracy.
   ---------
   The natural enclosure is extremely tight (radius ~ 2^{-prec}).  The result
   is returned at full precision; callers requiring a coarser bound may round
   after the fact.
*/

#include <flint/acb.h>
#include <flint/acb_poly.h>
#include <flint/acb_calc.h>
#include <flint/mag.h>
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Parameters threaded through the acb_calc / acb_calc_cauchy_bound callbacks. */
typedef struct {
  arb_t b;
  arb_t h;
  double mu;
} buthe_winf_params;

/*
   Direct factored evaluation of g(t) = E(t) * A(t) * B(t) on an arbitrary acb
   ball.  Used on [t0, T] and by the Cauchy bound on the circle |t| = R0, where
   the ball is bounded away from t = 0.  When order >= 1 the integrator requires
   analyticity on the ball; if a denominator (1 - e^{-2t}, cosh(h t/2), or t)
   straddles 0 the ball may contain a pole, so we return an indeterminate value
   and let the integrator subdivide.
*/
static int
buthe_winf_g(acb_ptr out, const acb_t inp, void *param, slong order, slong prec)
{
  buthe_winf_params *P = (buthe_winf_params *) param;
  acb_t t, E, A, B, tmp, arg, ch, sn, coef;
  arb_t pi;
  int pole;

  acb_init(t); acb_init(E); acb_init(A); acb_init(B);
  acb_init(tmp); acb_init(arg); acb_init(ch); acb_init(sn); acb_init(coef);
  arb_init(pi);
  arb_const_pi(pi, prec);
  acb_set(t, inp);

  /* E = exp(-(1/2 + mu) t) */
  arb_set_d(acb_realref(coef), 0.5 + P->mu);
  acb_mul(E, t, coef, prec);
  acb_neg(E, E);
  acb_exp(E, E, prec);

  /* A = 1 / (1 - e^{-2t}) */
  acb_mul_2exp_si(tmp, t, 1);
  acb_neg(tmp, tmp);
  acb_exp(tmp, tmp, prec);
  acb_sub_ui(tmp, tmp, 1, prec);
  acb_neg(tmp, tmp);                 /* 1 - e^{-2t} */
  pole = acb_contains_zero(tmp);
  acb_inv(A, tmp, prec);

  /* B = b/pi - sin(b t) / (pi t cosh(h t / 2)) */
  acb_mul_arb(arg, t, P->b, prec);
  acb_sin(sn, arg, prec);
  acb_mul_arb(arg, t, P->h, prec);
  acb_mul_2exp_si(arg, arg, -1);
  acb_cosh(ch, arg, prec);
  pole = pole || acb_contains_zero(ch) || acb_contains_zero(t);
  acb_mul(tmp, t, ch, prec);
  acb_div(tmp, sn, tmp, prec);
  acb_div_arb(tmp, tmp, pi, prec);
  acb_set_arb(B, P->b);
  acb_div_arb(B, B, pi, prec);
  acb_sub(B, B, tmp, prec);

  acb_mul(out, E, A, prec);
  acb_mul(out, out, B, prec);

  if (order >= 1 && pole)
    acb_indeterminate(out);

  arb_clear(pi);
  acb_clear(t); acb_clear(E); acb_clear(A); acb_clear(B);
  acb_clear(tmp); acb_clear(arg); acb_clear(ch); acb_clear(sn); acb_clear(coef);
  return 0;
}

/*
   Taylor series of g at t = 0, to length M, written into gc[0 .. M-1].
   g = E * Q / v with
       E(t) = e^{-(1/2+mu) t},
       v(t) = (1 - e^{-2t}) / (2t)            (entire, v(0) = 1),
       Q(t) = B(t) / (2t)                     (entire, Q(0) = 0),
       B(t) = (1/pi) ( b - (sin(b t)/t) / cosh(h t / 2) ).
   Every 1/t is realised as the exact shift-right of a series with a provably
   zero constant term, so no 0/0 arises and the coefficients are rigorous.
*/
static void
buthe_winf_series_at0(acb_poly_t gc, slong M, const arb_t b, const arb_t h,
                      double mu, slong prec)
{
  acb_poly_t X, E, v, sn, cc, Sc, Bb, Q, EQ;
  arb_t pi;
  acb_t c, c0;

  acb_poly_init(X); acb_poly_init(E); acb_poly_init(v);
  acb_poly_init(sn); acb_poly_init(cc); acb_poly_init(Sc);
  acb_poly_init(Bb); acb_poly_init(Q); acb_poly_init(EQ);
  arb_init(pi); arb_const_pi(pi, prec);
  acb_init(c); acb_init(c0);

  acb_poly_set_coeff_si(X, 1, 1);          /* X = t */

  /* E = exp(-(1/2 + mu) X) */
  arb_set_d(acb_realref(c), -(0.5 + mu));
  arb_zero(acb_imagref(c));
  acb_poly_scalar_mul(E, X, c, prec);
  acb_poly_exp_series(E, E, M, prec);

  /* v = (1 - e^{-2X}) / (2X).  Build -e^{-2X}, set the (provably 0) constant
     term of 1 - e^{-2X} to exact 0, shift right by 1, divide by 2. */
  acb_set_si(c, -2);
  acb_poly_scalar_mul(v, X, c, prec);
  acb_poly_exp_series(v, v, M + 1, prec);  /* e^{-2X} */
  acb_poly_neg(v, v);                      /* -e^{-2X}, constant term -1 */
  acb_zero(c);
  acb_poly_set_coeff_acb(v, 0, c);         /* 1 - e^{-2X}: constant term = 0 */
  acb_poly_shift_right(v, v, 1);           /* (1 - e^{-2X}) / X */
  acb_poly_truncate(v, M);
  acb_poly_scalar_mul_2exp_si(v, v, -1);   /* (1 - e^{-2X}) / (2X), v(0) = 1 */

  /* sin(b X) (constant term 0) and cosh(h X / 2) (constant term 1) */
  arb_set(acb_realref(c), b);
  arb_zero(acb_imagref(c));
  acb_poly_scalar_mul(sn, X, c, prec);
  acb_poly_sin_series(sn, sn, M + 1, prec);
  arb_set(acb_realref(c), h);
  arb_zero(acb_imagref(c));
  acb_poly_scalar_mul(cc, X, c, prec);
  acb_poly_scalar_mul_2exp_si(cc, cc, -1);
  acb_poly_cosh_series(cc, cc, M, prec);

  /* Sc = (sin(b X) / X) / cosh(h X / 2),  Sc(0) = b */
  acb_poly_shift_right(Sc, sn, 1);
  acb_poly_truncate(Sc, M);
  acb_poly_div_series(Sc, Sc, cc, M, prec);

  /* B = (1/pi) (b - Sc),  B(0) = 0 */
  acb_poly_neg(Bb, Sc);
  acb_poly_get_coeff_acb(c0, Bb, 0);
  arb_set(acb_realref(c), b);
  arb_zero(acb_imagref(c));
  acb_add(c0, c0, c, prec);
  acb_poly_set_coeff_acb(Bb, 0, c0);       /* b - Sc */
  arb_inv(acb_realref(c), pi, prec);
  arb_zero(acb_imagref(c));
  acb_poly_scalar_mul(Bb, Bb, c, prec);
  acb_zero(c);
  acb_poly_set_coeff_acb(Bb, 0, c);        /* B: constant term = 0 */

  /* Q = B / (2X) */
  acb_poly_shift_right(Q, Bb, 1);
  acb_poly_truncate(Q, M);
  acb_poly_scalar_mul_2exp_si(Q, Q, -1);

  /* g = E * Q / v */
  acb_poly_mullow(EQ, E, Q, M, prec);
  acb_poly_div_series(gc, EQ, v, M, prec);

  acb_clear(c); acb_clear(c0);
  arb_clear(pi);
  acb_poly_clear(X); acb_poly_clear(E); acb_poly_clear(v);
  acb_poly_clear(sn); acb_poly_clear(cc); acb_poly_clear(Sc);
  acb_poly_clear(Bb); acb_poly_clear(Q); acb_poly_clear(EQ);
}

void
buthe_winf_integral(arb_t res, const arb_t b, const arb_t h, double mu, slong prec)
{
  buthe_winf_params P;
  arb_t t0, R0, alpha, T, total, tmp, pi;
  slong wprec, M;

  /* Use a little extra internal precision to absorb rounding in intermediate steps. */
  wprec = prec + 32;

  arb_init(P.b); arb_init(P.h);
  arb_set(P.b, b); arb_set(P.h, h);
  P.mu = mu;

  arb_init(t0); arb_init(R0); arb_init(alpha);
  arb_init(T); arb_init(total); arb_init(tmp); arb_init(pi);
  arb_const_pi(pi, wprec);
  arb_zero(total);

  /* t0 = 1/32.  R0 = 0.9 * (pi/h): the nearest singularity of 1/cosh(ht/2)
     to the origin is the pole at t = i*pi/h, so any R0 < pi/h keeps the disc
     inside the analyticity strip.  Series terms past M are bounded by a Cauchy
     estimate with ratio t0/R0 = h/(28.8*pi) ~ h/90.5.  That ratio is small (a
     tight enclosure) for the h actually used: the h < 2*pi*b/5 guard caps
     h <= 8 for degree <= 9, giving ratio ~ 0.09.  As h approaches ~90 the ratio
     approaches 1 and the Cauchy tail degrades to a fail-safe non-finite ball
     (which can never false-confirm); h beyond ~90 is out of the supported
     range. */
  // R0 < pi/h keeps the disc inside the analyticity strip of 1/cosh(h t/2)
  arb_set_ui(t0, 1);
  arb_div_ui(t0, t0, 32, wprec);
  arb_const_pi(R0, wprec);
  arb_div(R0, R0, h, wprec);
  arb_mul_ui(R0, R0, 9, wprec);
  arb_div_ui(R0, R0, 10, wprec);
  arb_set_d(alpha, 0.5 + mu);

  /* Enough terms that the Cauchy tail (t0/R0)^M is far below 2^{-prec}. */
  M = prec / 3 + 48;
  if (M < 64)
    M = 64;

  /* ---------------- [0, t0] : Taylor series at 0 ---------------- */
  {
    acb_poly_t gc;
    arb_t acc, pw, gk, Cg, ratio, rem, one;
    acb_t center, ck;
    mag_t mrem;
    slong k;

    acb_poly_init(gc);
    buthe_winf_series_at0(gc, M, b, h, mu, wprec);

    arb_init(acc); arb_init(pw); arb_init(gk);
    arb_init(Cg); arb_init(ratio); arb_init(rem); arb_init(one);
    acb_init(center); acb_init(ck);
    mag_init(mrem);

    /* integrate the truncated series:  sum_{k=0}^{M-1} g_k t0^{k+1}/(k+1) */
    arb_zero(acc);
    arb_set(pw, t0);                       /* t0^{k+1}, starting k = 0 */
    for (k = 0; k < M; k++) {
      acb_poly_get_coeff_acb(ck, gc, k);
      arb_set(gk, acb_realref(ck));        /* g is real on the real axis */
      arb_mul(tmp, gk, pw, wprec);
      arb_div_ui(tmp, tmp, (ulong) (k + 1), wprec);
      arb_add(acc, acc, tmp, wprec);
      arb_mul(pw, pw, t0, wprec);
    }

    /* Cauchy tail bound:  |remainder| <= t0 * Cg * (t0/R0)^M / (1 - t0/R0),
       where Cg is a rigorous bound on |g| over the disc |t| <= R0. */
    acb_zero(center);
    acb_calc_cauchy_bound(Cg, buthe_winf_g, &P, center, R0, 16, wprec);
    arb_div(ratio, t0, R0, wprec);
    arb_pow_ui(rem, ratio, (ulong) M, wprec);
    arb_mul(rem, rem, Cg, wprec);
    arb_mul(rem, rem, t0, wprec);
    arb_set_ui(one, 1);
    arb_sub(one, one, ratio, wprec);
    arb_div(rem, rem, one, wprec);
    arb_get_mag(mrem, rem);
    arb_add_error_mag(acc, mrem);

    arb_add(total, total, acc, wprec);

    mag_clear(mrem);
    acb_clear(center); acb_clear(ck);
    arb_clear(acc); arb_clear(pw); arb_clear(gk);
    arb_clear(Cg); arb_clear(ratio); arb_clear(rem); arb_clear(one);
    acb_poly_clear(gc);
  }

  /* ---------------- choose T ---------------- */
  /* T ~ (prec + 40) ln2 / alpha + 12 makes the closed-form tail below
     2^{-prec}; correctness does not depend on this estimate being tight. */
  {
    arb_t ln2;
    arb_init(ln2);
    arb_const_log2(ln2, wprec);
    arb_mul_ui(T, ln2, (ulong) (prec + 40), wprec);
    arb_div(T, T, alpha, wprec);
    arb_add_ui(T, T, 12, wprec);
    arb_clear(ln2);
  }

  /* ---------------- [t0, T] : rigorous integrator ---------------- */
  {
    acb_t a, bb, ires;
    mag_t tol;
    acb_calc_integrate_opt_t opt;

    acb_init(a); acb_init(bb); acb_init(ires);
    mag_init(tol);
    arb_set(acb_realref(a), t0);
    arb_set(acb_realref(bb), T);
    mag_set_ui_2exp_si(tol, 1, -(prec + 8));
    acb_calc_integrate_opt_init(opt);
    /* status != 0 means the integrator did not converge to tol; the result is
       still a rigorous enclosure (FLINT contract), just wider than requested. */
    int status = acb_calc_integrate(ires, buthe_winf_g, &P, a, bb, prec + 8, tol, opt, wprec);
    (void) status;
    arb_add(total, total, acb_realref(ires), wprec);

    mag_clear(tol);
    acb_clear(a); acb_clear(bb); acb_clear(ires);
  }

  /* ---------------- [T, infty) : closed-form bound ---------------- */
  {
    arb_t bnd, denom, expo, brk, zero;
    mag_t mb;

    arb_init(bnd); arb_init(denom); arb_init(expo); arb_init(brk); arb_init(zero);
    mag_init(mb);

    /* 1 - e^{-2T} */
    arb_mul_2exp_si(tmp, T, 1);
    arb_neg(tmp, tmp);
    arb_exp(tmp, tmp, wprec);
    arb_set_ui(denom, 1);
    arb_sub(denom, denom, tmp, wprec);
    /* e^{-alpha T} / alpha */
    arb_mul(tmp, alpha, T, wprec);
    arb_neg(tmp, tmp);
    arb_exp(expo, tmp, wprec);
    arb_div(expo, expo, alpha, wprec);
    /* b/pi + 1/(pi T) */
    arb_div(brk, b, pi, wprec);
    arb_mul(tmp, pi, T, wprec);
    arb_inv(tmp, tmp, wprec);
    arb_add(brk, brk, tmp, wprec);
    /* bnd = (e^{-alpha T}/alpha) / (1 - e^{-2T}) * (b/pi + 1/(pi T)) */
    arb_div(bnd, expo, denom, wprec);
    arb_mul(bnd, bnd, brk, wprec);

    arb_get_mag(mb, bnd);
    arb_zero(zero);
    arb_add_error_mag(zero, mb);
    arb_add(total, total, zero, wprec);

    mag_clear(mb);
    arb_clear(bnd); arb_clear(denom); arb_clear(expo); arb_clear(brk); arb_clear(zero);
  }

  /* I = 2 * (sum of the three pieces) */
  arb_mul_2exp_si(total, total, 1);

  arb_set(res, total);

  arb_clear(t0); arb_clear(R0); arb_clear(alpha);
  arb_clear(T); arb_clear(total); arb_clear(tmp); arb_clear(pi);
  arb_clear(P.b); arb_clear(P.h);
}

#ifdef __cplusplus
}
#endif
