#include <assert.h>
#include <inttypes.h>
#include <flint/acb_poly.h>
#include <flint/fmpz_poly.h>
#include "glfunc.h"
#include "sym_power.h"

static long g_ainvs[5];
static long g_bad;
static long g_apbad;
static int g_sym;

static long mod_p(long value, long p)
{
  return ((value % p) + p) % p;
}

static long ap_good(long p)
{
  long a1 = g_ainvs[0], a2 = g_ainvs[1], a3 = g_ainvs[2];
  long a4 = g_ainvs[3], a6 = g_ainvs[4];
  long count = 1;
  for (long x = 0; x < p; x++) {
    long x2 = mod_p(x * x, p), rhs = mod_p(x2 * x, p);
    rhs = mod_p(rhs + a2 * x2, p);
    rhs = mod_p(rhs + a4 * x, p);
    rhs = mod_p(rhs + a6, p);
    for (long y = 0; y < p; y++)
      if (mod_p(y * y + a1 * x * y + a3 * y, p) == rhs)
        count++;
  }
  return p + 1 - count;
}

static void set_sym2_factor(acb_poly_t poly, long ap, uint64_t p)
{
  fmpz_poly_t factor;
  fmpz_t coefficient;
  acb_t value;
  fmpz_poly_init(factor);
  fmpz_init(coefficient);
  acb_init(value);
  sym_power_lpoly(factor, ap, p, 2);
  for (slong i = 0; i <= fmpz_poly_degree(factor); i++) {
    fmpz_poly_get_coeff_fmpz(coefficient, factor, i);
    acb_set_fmpz(value, coefficient);
    acb_poly_set_coeff_acb(poly, i, value);
  }
  acb_clear(value);
  fmpz_clear(coefficient);
  fmpz_poly_clear(factor);
}

static void local_factor(acb_poly_t poly, uint64_t p, int degree,
                         int64_t prec, void *param)
{
  (void)degree;
  (void)prec;
  (void)param;
  acb_poly_zero(poly);
  long prime = (long)p;
  if (prime == g_bad) {
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, g_sym == 2 ? -1 : -g_apbad);
    return;
  }

  long ap = ap_good(prime);
  if (g_sym == 2) {
    set_sym2_factor(poly, ap, p);
  } else {
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, -ap);
    acb_poly_set_coeff_si(poly, 2, prime);
  }
}

static Lerror_t run(long ainvs[5], uint64_t conductor, long bad, long apbad,
                    int sym, uint64_t degree, double normalisation,
                    double *mus, int supplied_rank, acb_t sign_out,
                    int64_t *rank_out)
{
  for (int i = 0; i < 5; i++)
    g_ainvs[i] = ainvs[i];
  g_bad = bad;
  g_apbad = apbad;
  g_sym = sym;

  Lerror_t ecode = ERR_SUCCESS;
  Lparams_t params = {0};
  params.degree = degree;
  params.conductor = conductor;
  params.normalisation = normalisation;
  params.mus = mus;
  params.target_prec = DEFAULT_TARGET_PREC;
  params.self_dual = YES;
  params.rank = supplied_rank;
  params.cache_dir = NULL;

  Lfunc_t L = Lfunc_init_advanced(&params, &ecode);
  ecode |= Lfunc_use_all_lpolys(L, local_factor, NULL);
  ecode |= Lfunc_compute(L);
  acb_set(sign_out, Lfunc_sign(L));
  *rank_out = Lfunc_rank(L);
  Lfunc_clear(L);
  return ecode;
}

int main(void)
{
  long a11[5] = {0, -1, 1, -10, -20};
  long a37[5] = {0, 0, 1, -1, 0};
  double mus2[2] = {0, 1};
  double mus3[3] = {0, 1, 0};
  acb_t sign;
  int64_t rank;
  acb_init(sign);

  Lerror_t seed = run(a11, 121, 11, 1, 2, 3, 1.0, mus3, 0,
                      sign, &rank);
  assert(!fatal_error(seed));
  assert(rank == 0);
  assert((seed & (ERR_INSUFF_EULER | ERR_NO_RANK | ERR_CONFLICT_RANK |
                  ERR_SOME_DATA | ERR_ZERO_PREC | ERR_DBL_ZERO)) == 0);

  Lerror_t after = run(a37, 37, 37, -1, 1, 2, 0.5, mus2, DK,
                       sign, &rank);
  assert(rank == 1);
  assert(arb_is_negative(acb_realref(sign)));
  assert(arb_contains_zero(acb_imagref(sign)));
  assert(after == ERR_SUCCESS);

  Lerror_t conflict = run(a11, 121, 11, 1, 2, 3, 1.0, mus3, 1,
                          sign, &rank);
  assert(!fatal_error(conflict));
  assert(rank == 1);
  assert((conflict & (ERR_CONFLICT_RANK | ERR_RH_ERROR)) ==
         (ERR_CONFLICT_RANK | ERR_RH_ERROR));
  assert((conflict & (ERR_INSUFF_EULER | ERR_NO_RANK | ERR_SOME_DATA |
                      ERR_ZERO_PREC | ERR_DBL_ZERO)) == 0);

  acb_clear(sign);
  return 0;
}
