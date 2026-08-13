#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <flint/fmpq.h>
#include "glfunc_internals.h"

static void set_rational(arb_t x, slong numerator, ulong denominator)
{
  fmpq_t q;
  fmpq_init(q);
  fmpq_set_si(q, numerator, denominator);
  arb_set_fmpq(x, q, DEFAULT_TARGET_PREC);
  fmpq_clear(q);
}

static void set_interval(arb_t x, slong lo_num, ulong lo_den,
                         slong hi_num, ulong hi_den)
{
  arb_t lo, hi;
  arb_init(lo);
  arb_init(hi);
  set_rational(lo, lo_num, lo_den);
  set_rational(hi, hi_num, hi_den);
  arb_union(x, lo, hi, DEFAULT_TARGET_PREC);
  arb_clear(lo);
  arb_clear(hi);
}

static void set_sign(acb_t sign, slong re_lo_num, ulong re_lo_den,
                     slong re_hi_num, ulong re_hi_den,
                     slong im_lo_num, ulong im_lo_den,
                     slong im_hi_num, ulong im_hi_den)
{
  set_interval(acb_realref(sign), re_lo_num, re_lo_den,
               re_hi_num, re_hi_den);
  set_interval(acb_imagref(sign), im_lo_num, im_lo_den,
               im_hi_num, im_hi_den);
}

static void assert_evidence(turing_count_evidence_t evidence,
                            uint64_t central_zeros, bool certified)
{
  assert(evidence.central_zeros == central_zeros);
  assert(evidence.certified == certified);
}

static void test_numeric_classifier(void)
{
  arb_t tcount;
  acb_t sign;
  arb_init(tcount);
  acb_init(sign);

  turing_count_evidence_t certified = {0, true};
  set_interval(tcount, 367722426057, 10000000000,
               393350136280, 10000000000);
  assert(turing_count_status_with_evidence(
             certified, tcount, 38, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_TOO_FEW);
  assert(turing_count_status_with_evidence(
             certified, tcount, 39, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_CONFIRMED);

  set_sign(sign, -3, 2, -1, 2, 0, 1, 0, 1);
  turing_count_evidence_t conflict = turing_make_count_evidence(
      1, ERR_CONFLICT_RANK, YES, sign, ERR_SUCCESS, ERR_NO_DATA);
  assert(turing_count_status_with_evidence(
             conflict, tcount, 39, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_UNCERTIFIED);

  arb_set_ui(tcount, 4);
  assert(turing_count_status_with_evidence(
             certified, tcount, 5, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_TOO_MANY);
  arb_set_ui(tcount, 5);
  assert(turing_count_status_with_evidence(
             certified, tcount, 5, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_CONFIRMED);
  arb_set_ui(tcount, 6);
  assert(turing_count_status_with_evidence(
             certified, tcount, 5, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_TOO_FEW);

  set_interval(tcount, 4, 1, 6, 1);
  assert(turing_count_status_with_evidence(
             conflict, tcount, 5, DEFAULT_TARGET_PREC) ==
         TURING_COUNT_UNCERTIFIED);

  acb_clear(sign);
  arb_clear(tcount);
}

static void test_constructor_signs(void)
{
  acb_t sign;
  acb_init(sign);

  set_sign(sign, -3, 2, -1, 2, 0, 1, 0, 1);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  1, true);

  set_sign(sign, -3, 2, -1, 2, -1, 10, 1, 10);
  assert(!acb_is_real(sign));
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  1, true);

  set_sign(sign, -3, 2, 1, 2, 0, 1, 0, 1);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);

  set_sign(sign, 1, 2, 3, 2, 0, 1, 0, 1);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);

  set_sign(sign, -3, 2, -1, 2, 1, 10, 2, 10);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);

  set_sign(sign, -3, 2, -1, 2, 0, 1, 0, 1);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, NO, sign,
                      ERR_SUCCESS, ERR_SUCCESS),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_SUCCESS, DK, sign,
                      ERR_SUCCESS, ERR_SUCCESS),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);

  set_sign(sign, 1, 2, 3, 2, 0, 1, 0, 1);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, true);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, NO, sign,
                      ERR_SUCCESS, ERR_SUCCESS),
                  0, true);

  acb_clear(sign);
}

static void test_constructor_blockers(void)
{
  static const Lerror_t warnings[] = {
      ERR_SOME_DATA, ERR_ZERO_PREC, ERR_RH_ERROR, ERR_INSUFF_EULER,
      ERR_NO_RANK, ERR_CONFLICT_RANK, ERR_DBL_ZERO, ERR_SPEC_PREC,
      ERR_G_OUTFILE};
  acb_t sign;
  acb_init(sign);
  set_sign(sign, 1, 2, 3, 2, 0, 1, 0, 1);

  assert_evidence(turing_make_count_evidence(
                      0, ERR_CONFLICT_RANK, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      1, ERR_CONFLICT_RANK, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_NO_RANK, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      2, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      DK, ERR_SUCCESS, YES, sign,
                      ERR_SUCCESS, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, 2, sign,
                      ERR_SUCCESS, ERR_SUCCESS),
                  0, false);

  for (size_t i = 0; i < sizeof(warnings) / sizeof(warnings[0]); i++) {
    assert_evidence(turing_make_count_evidence(
                        0, ERR_SUCCESS, YES, sign,
                        warnings[i], ERR_NO_DATA),
                    0, false);
    assert_evidence(turing_make_count_evidence(
                        0, ERR_SUCCESS, NO, sign,
                        warnings[i], ERR_SUCCESS),
                    0, false);
    assert_evidence(turing_make_count_evidence(
                        0, ERR_SUCCESS, NO, sign,
                        ERR_SUCCESS, warnings[i]),
                    0, false);
    assert_evidence(turing_make_count_evidence(
                        0, ERR_SUCCESS, DK, sign,
                        ERR_SUCCESS, warnings[i]),
                    0, false);
    assert_evidence(turing_make_count_evidence(
                        0, ERR_SUCCESS, YES, sign,
                        ERR_SUCCESS, warnings[i]),
                    0, true);
  }

  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, YES, sign,
                      ERR_UPSAMPLE, ERR_NO_DATA),
                  0, false);
  assert_evidence(turing_make_count_evidence(
                      0, ERR_SUCCESS, NO, sign,
                      ERR_SUCCESS, ERR_STAT_POINT),
                  0, false);

  acb_clear(sign);
}

static void test_uncertified_early_gate(void)
{
  Lfunc L = {0};
  arb_t sentinel;
  arb_init(L.imint);
  arb_init(sentinel);
  arb_set_ui(sentinel, 7);
  arb_set(L.imint, sentinel);
  L.zeros[0] = NULL;
  L.zeros[1] = (arb_t *)(uintptr_t)1;

  turing_count_evidence_t uncertified = {0, false};
  assert(turing_check_RH(&L, uncertified, DEFAULT_TARGET_PREC) ==
         ERR_RH_ERROR);
  assert(arb_equal(L.imint, sentinel));

  arb_clear(sentinel);
  arb_clear(L.imint);
}

int main(int argc, char **argv)
{
  if (argc == 1 || strcmp(argv[1], "numeric") == 0)
    test_numeric_classifier();
  if (argc == 1 || strcmp(argv[1], "signs") == 0)
    test_constructor_signs();
  if (argc == 1 || strcmp(argv[1], "blockers") == 0)
    test_constructor_blockers();
  if (argc == 1 || strcmp(argv[1], "gate") == 0)
    test_uncertified_early_gate();
  return 0;
}
