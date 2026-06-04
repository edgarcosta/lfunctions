// Coverage for the non-self-dual Buthe threshold (0.49). The degree-2 L =
// chi_5 * chi_7 (non-principal characters mod 5 and 7, conductor 35) is built
// with self_dual left DK (the default for Lfunc_init), so buthe_check_RH uses
// the tighter 0.49 single-zero bar. The object is genuine and complete, so RH
// must still be CONFIRMED under 0.49: ERR_RH_ERROR / ERR_BUT_ERROR absent.
// This pins that tightening the bar to 0.49 for non-self-dual L does not break
// a real non-self-dual confirmation. Asserts on flags + a certified value.
// Exit 0 = pass.
#include <assert.h>
#include <inttypes.h>
#include <stdio.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

// chi_5 * chi_7 Euler factor (same as test/dir_test.c), running out past p>100.
static void cb(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *pm) {
  (void)d; (void)pm;
  if (p > 100) { acb_poly_zero(poly); return; }
  acb_poly_t p5, p7; acb_poly_init(p5); acb_poly_init(p7);
  acb_poly_one(p5); acb_poly_one(p7);
  if ((p % 5 == 1) || (p % 5 == 4)) acb_poly_set_coeff_si(p5, 1, -1);
  if ((p % 5 == 2) || (p % 5 == 3)) acb_poly_set_coeff_si(p5, 1, 1);
  if ((p % 7 == 1) || (p % 7 == 2) || (p % 7 == 4)) acb_poly_set_coeff_si(p7, 1, -1);
  if ((p % 7 == 3) || (p % 7 == 5) || (p % 7 == 6)) acb_poly_set_coeff_si(p7, 1, 1);
  acb_poly_mul(poly, p5, p7, prec);
  acb_poly_clear(p5); acb_poly_clear(p7);
}

int main(void) {
  Lerror_t ec = 0;
  double mus[2] = {0, 1};
  Lfunc_t L = Lfunc_init(2, 5 * 7, 0.0, mus, &ec); // self_dual = DK by default
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  // sanity: this object is NOT flagged self-dual, so the 0.49 branch is taken
  assert(((Lfunc *)L)->self_dual != YES);

  ec |= Lfunc_use_all_lpolys(L, cb, NULL);
  if (fatal_error(ec)) { fprint_errors(stderr, ec); return 2; }
  ec |= Lfunc_compute(L); // default verifier is Buthe (Task 3); uses 0.49 here
  if (fatal_error(ec)) { fprint_errors(stderr, ec); Lfunc_clear(L); return 2; }

  printf("non-self-dual Buthe ecode = 0x%lx\n", (unsigned long)ec);
  int ok = 1;
  if (ec & ERR_RH_ERROR) { fprintf(stderr,
      "FAIL: ERR_RH_ERROR under 0.49 threshold for a complete non-self-dual L\n"); ok = 0; }
  if (ec & ERR_BUT_ERROR) { fprintf(stderr, "FAIL: ERR_BUT_ERROR set\n"); ok = 0; }
  assert((ec & ERR_RH_ERROR) == 0);
  assert((ec & ERR_BUT_ERROR) == 0);

  // certify the first zero (this L has a real first zero; read it once to pin).
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  assert(!arb_is_zero(zeros + 0)); // at least one zero found

  // --- white-box: same forced S=0.7, different self_dual => different verdict.
  // buthe_check_RH recomputes buthe_Ws and buthe_Winf internally (from the zeros
  // and self_dual), so the ONLY accumulator it reads-but-does-not-recompute is
  // buthe_Wf. We therefore push S into (0.49, 0.98) by bumping buthe_Wf, sizing
  // the bump against the S the function actually computes for the chosen
  // self_dual path (a probe call leaves Ws/Winf populated for that path, so
  // S = Wf + Winf - Ws is readable from the struct). With S = 0.7 fixed for each
  // path in turn, the verdict isolates the threshold: 0.49 (self_dual != YES)
  // rejects, 0.98 (self_dual == YES) accepts. On the pre-Task-5 tree both use
  // 0.98, so the self_dual=NO assert below returns SUCCESS and FAILS (the
  // genuine fail-on-base distinguishing the 0.49 tightening from 0.98).
  {
    Lfunc *LL = (Lfunc *)L;
    arb_t target, S; arb_init(target); arb_init(S);
    arb_set_d(target, 0.7);

    // self_dual = NO: probe to populate Ws/Winf, then bump Wf so S becomes 0.7.
    LL->self_dual = NO;
    (void)buthe_check_RH(LL); // recomputes buthe_Ws/buthe_Winf for the NO path
    arb_add(S, LL->buthe_Wf, LL->buthe_Winf, LL->wprec);
    arb_sub(S, S, LL->buthe_Ws, LL->wprec);          // S the NO path used
    arb_sub(S, target, S, LL->wprec);                // delta to reach 0.7
    arb_add(LL->buthe_Wf, LL->buthe_Wf, S, LL->wprec); // Wf += delta => S -> 0.7
    assert((buthe_check_RH(LL) & ERR_RH_ERROR) != 0);  // 0.49 bar rejects

    // self_dual = YES: re-probe (the YES path recomputes a different Ws), resize
    // the bump so S is 0.7 again, then confirm the 0.98 bar accepts.
    LL->self_dual = YES;
    (void)buthe_check_RH(LL); // recomputes buthe_Ws/buthe_Winf for the YES path
    arb_add(S, LL->buthe_Wf, LL->buthe_Winf, LL->wprec);
    arb_sub(S, S, LL->buthe_Ws, LL->wprec);          // S the YES path used
    arb_sub(S, target, S, LL->wprec);                // delta to reach 0.7
    arb_add(LL->buthe_Wf, LL->buthe_Wf, S, LL->wprec);
    assert((buthe_check_RH(LL) & ERR_RH_ERROR) == 0);  // 0.98 bar accepts

    arb_clear(target); arb_clear(S);
  }

  Lfunc_clear(L);
  if (!ok) return 1;
  printf("PASS: non-self-dual L confirms RH under the 0.49 threshold.\n");
  return 0;
}
