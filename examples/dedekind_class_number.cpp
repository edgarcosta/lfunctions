// Self-test examples for number-field quotients:
// M_K(s) = zeta_K(s) / zeta(s), checked at s = 1 against the analytic
// class number formula and against zeros.
//
// The engine does not support the pole of zeta_K directly, so this example
// verifies the Dedekind residue through the entire quotient M_K.
// PARI's lfuncreate(polynomial) computes zeta_K, so zero goldens for M_K are
// formed by omitting the Riemann-zeta zeros from the PARI list when they occur.
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/fmpz_poly.h>
#include <flint/mag.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

#include "glfunc.h"

using std::map;
using std::uint64_t;
using std::vector;

static const size_t MAX_TEST_DEGREE = 8;
static const size_t MAX_TEST_ZEROS = 9;

struct nf_ctx {
  vector<long> defining_poly;            // ascending coefficients
  map<uint64_t, vector<long> > bad_lpolys; // ascending M_K local polynomial coeffs
};

struct bad_factor {
  uint64_t p;
  const long *coeffs;
  size_t ncoeffs;
};

struct dedekind_testcase {
  const char *label;
  const long *defining_poly;
  size_t ncoeffs;
  const bad_factor *bad_factors;
  size_t nbad_factors;
  int degree;                    // degree of M_K = zeta_K / zeta
  uint64_t conductor;
  const double mus[MAX_TEST_DEGREE];  // gamma shifts for M_K
  int r1, r2;
  uint64_t class_number;
  const char *regulator;
  uint64_t torsion_order;
  uint64_t disc_abs;
  int64_t rank;
  const char *epsilon_re, *epsilon_im;
  const char *zeros0[MAX_TEST_ZEROS];
};

static ulong coeff_mod(long c, ulong p)
{
  if (c >= 0)
    return ulong(c) % p;
  const ulong a = ulong(-c) % p;
  return a == 0 ? 0 : p - a;
}

static void fail_missing_bad_factor(uint64_t p)
{
  fprintf(stderr, "dedekind_class_number: missing bad local factor for prime %llu\n",
          (unsigned long long) p);
  std::abort();
}

static void factor_degrees_modp(vector<slong>& degrees,
                                const vector<long>& polynomial,
                                uint64_t p)
{
  nmod_poly_t pol;
  nmod_poly_init2(pol, ulong(p), slong(polynomial.size()));
  for (size_t i = 0; i < polynomial.size(); ++i)
    nmod_poly_set_coeff_ui(pol, slong(i), coeff_mod(polynomial[i], ulong(p)));

  nmod_poly_factor_t fac;
  nmod_poly_factor_init(fac);
  nmod_poly_factor(fac, pol);

  degrees.clear();
  for (slong i = 0; i < fac->num; ++i) {
    if (fac->exp[i] != 1) {
      nmod_poly_factor_clear(fac);
      nmod_poly_clear(pol);
      fail_missing_bad_factor(p);
    }
    degrees.push_back(nmod_poly_degree(fac->p + i));
  }
  std::sort(degrees.begin(), degrees.end());

  nmod_poly_factor_clear(fac);
  nmod_poly_clear(pol);
}

static void mk_lpoly_from_degrees(fmpz_poly_t out, const vector<slong>& degrees)
{
  fmpz_poly_t prod, factor, divisor, quotient;
  fmpz_poly_init(prod);
  fmpz_poly_init(factor);
  fmpz_poly_init(divisor);
  fmpz_poly_init(quotient);

  fmpz_poly_one(prod);
  for (slong d : degrees) {
    fmpz_poly_zero(factor);
    fmpz_poly_set_coeff_si(factor, 0, 1);
    fmpz_poly_set_coeff_si(factor, d, -1);
    fmpz_poly_mul(prod, prod, factor);
  }

  fmpz_poly_set_coeff_si(divisor, 0, 1);
  fmpz_poly_set_coeff_si(divisor, 1, -1);
  fmpz_poly_divexact(quotient, prod, divisor);
  fmpz_poly_set(out, quotient);

  fmpz_poly_clear(prod);
  fmpz_poly_clear(factor);
  fmpz_poly_clear(divisor);
  fmpz_poly_clear(quotient);
}

static void assert_lpoly_coeffs(const vector<slong>& degrees, const vector<long>& expected)
{
  fmpz_poly_t pol;
  fmpz_poly_init(pol);
  mk_lpoly_from_degrees(pol, degrees);

  assert(fmpz_poly_length(pol) == slong(expected.size()));
  for (size_t i = 0; i < expected.size(); ++i)
    assert(fmpz_poly_get_coeff_si(pol, slong(i)) == expected[i]);

  fmpz_poly_clear(pol);
}

static void assert_lpoly_sanity_checks()
{
  assert_lpoly_coeffs(vector<slong>{1, 1, 1}, vector<long>{1, -2, 1});
  assert_lpoly_coeffs(vector<slong>{1, 2}, vector<long>{1, 0, -1});
  assert_lpoly_coeffs(vector<slong>{3}, vector<long>{1, 1, 1});
}

static void set_poly_from_coeffs(acb_poly_t poly, const vector<long>& coeffs)
{
  acb_poly_zero(poly);
  for (size_t i = 0; i < coeffs.size(); ++i)
    acb_poly_set_coeff_si(poly, slong(i), coeffs[i]);
}

static void mk_lpoly_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  (void) d;
  nf_ctx *ctx = static_cast<nf_ctx *>(param);
  const auto bad = ctx->bad_lpolys.find(p);
  if (bad != ctx->bad_lpolys.end()) {
    set_poly_from_coeffs(poly, bad->second);
    return;
  }

  vector<slong> degrees;
  factor_degrees_modp(degrees, ctx->defining_poly, p);

  fmpz_poly_t f;
  fmpz_poly_init(f);
  mk_lpoly_from_degrees(f, degrees);
  acb_poly_set_fmpz_poly(poly, f, prec);
  fmpz_poly_clear(f);
}

static void analytic_class_number_rhs(acb_t rhs, const dedekind_testcase& c,
                                      slong prec)
{
  arb_t regulator, pi, sqrt_disc, real, term;
  arb_init(regulator);
  arb_init(pi);
  arb_init(sqrt_disc);
  arb_init(real);
  arb_init(term);

  arb_set_str(regulator, c.regulator, prec);
  arb_const_pi(pi, prec);
  arb_sqrt_ui(sqrt_disc, c.disc_abs, prec);

  // 2^r1 (2*pi)^r2 h R / (w sqrt(|D|)).
  arb_set_ui(real, 1);
  arb_mul_2exp_si(real, real, c.r1);
  arb_mul_ui(term, pi, 2, prec);
  for (int i = 0; i < c.r2; ++i)
    arb_mul(real, real, term, prec);
  arb_mul_ui(real, real, c.class_number, prec);
  arb_mul(real, real, regulator, prec);
  arb_div_ui(real, real, c.torsion_order, prec);
  arb_div(real, real, sqrt_disc, prec);

  acb_zero(rhs);
  arb_set(acb_realref(rhs), real);

  arb_clear(regulator);
  arb_clear(pi);
  arb_clear(sqrt_disc);
  arb_clear(real);
  arb_clear(term);
}

static void assert_radius_below(arb_srcptr x, double tol)
{
  mag_t mtol;
  mag_init(mtol);
  mag_set_d(mtol, tol);
  assert(mag_cmp(arb_radref(x), mtol) < 0);
  mag_clear(mtol);
}

static void assert_arb_overlaps_str(arb_srcptr value, const char *s)
{
  arb_t ref;
  arb_init(ref);
  arb_set_str(ref, s, 300);
  arb_add_error_2exp_si(ref, -50);
  assert(arb_overlaps(value, ref));
  arb_clear(ref);
}

static void assert_acb_overlaps_str(acb_srcptr value, const char *re, const char *im)
{
  acb_t ref;
  acb_init(ref);
  arb_set_str(acb_realref(ref), re, 300);
  arb_set_str(acb_imagref(ref), im, 300);
  arb_add_error_2exp_si(acb_realref(ref), -50);
  arb_add_error_2exp_si(acb_imagref(ref), -50);
  assert(acb_overlaps(value, ref));
  acb_clear(ref);
}

static void assert_sign(Lfunc_t L, const dedekind_testcase& c)
{
  const slong prec = 300;
  arb_t abs, one;
  arb_init(abs);
  arb_init(one);
  acb_abs(abs, Lfunc_sign(L), prec);
  arb_one(one);
  assert(arb_overlaps(abs, one));
  arb_clear(abs);
  arb_clear(one);

  assert_acb_overlaps_str(Lfunc_sign(L), c.epsilon_re, c.epsilon_im);
}

static void assert_zeros(Lfunc_t L, const dedekind_testcase& c)
{
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  for (size_t i = 0; i < MAX_TEST_ZEROS && c.zeros0[i] != NULL; ++i) {
    assert(!arb_contains_zero(zeros + i));
    assert_arb_overlaps_str(zeros + i, c.zeros0[i]);
  }
}

static void print_zeros(Lfunc_t L, const dedekind_testcase& c)
{
  if (getenv("DEDEKIND_PRINT_ZEROS") == NULL)
    return;

  arb_srcptr zeros = Lfunc_zeros(L, 0);
  printf("dedekind_class_number: zeros %s", c.label);
  for (size_t i = 0; i < MAX_TEST_ZEROS && !arb_is_zero(zeros + i); ++i) {
    printf(" ");
    arb_printd(zeros + i, 30);
  }
  printf("\n");
  fflush(stdout);
}

static Lfunc_t init_case(const dedekind_testcase& c, nf_ctx& ctx, Lerror_t& ecode)
{
  ctx.defining_poly.assign(c.defining_poly, c.defining_poly + c.ncoeffs);
  ctx.bad_lpolys.clear();
  for (size_t i = 0; i < c.nbad_factors; ++i) {
    const bad_factor& bf = c.bad_factors[i];
    ctx.bad_lpolys[bf.p] = vector<long>(bf.coeffs, bf.coeffs + bf.ncoeffs);
  }

  static char cache_dir[] = ".";
  Lparams_t Lp;
  Lp.degree = c.degree;
  Lp.conductor = c.conductor;
  Lp.normalisation = 0.0;
  Lp.mus = const_cast<double *>(c.mus);
  Lp.target_prec = 160;
  Lp.wprec = 0;
  Lp.gprec = 0;
  Lp.self_dual = YES;
  Lp.rank = DK;
  Lp.cache_dir = cache_dir;

  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  assert(!fatal_error(ecode));

  ecode |= Lfunc_use_all_lpolys(L, mk_lpoly_callback, &ctx);
  assert(!fatal_error(ecode));

  ecode |= Lfunc_compute(L);
  assert(!fatal_error(ecode));
  assert(Lfunc_rank(L) == c.rank);
  return L;
}

static void run_case(const dedekind_testcase& c)
{
  printf("dedekind_class_number: running %s\n", c.label);
  fflush(stdout);
  nf_ctx ctx;
  Lerror_t ecode = 0;
  Lfunc_t L = init_case(c, ctx, ecode);
  assert(Lfunc_rank(L) == c.rank);
  assert_sign(L, c);

  acb_t MK1, rhs;
  acb_init(MK1);
  acb_init(rhs);
  analytic_class_number_rhs(rhs, c, 300);

  const Lerror_t sv_error = Lfunc_special_value(MK1, L, 1.0, 0.0);
  ecode |= sv_error;
  assert(!fatal_error(sv_error));
  assert(!fatal_error(ecode));
  assert((sv_error & ERR_SPEC_PREC) == 0);

  if (!acb_overlaps(MK1, rhs)) {
    printf("  M_K(1) = "); acb_printd(MK1, 40); printf("\n");
    printf("  RHS    = "); acb_printd(rhs, 40); printf("\n");
    fflush(stdout);
  }
  assert(acb_overlaps(MK1, rhs));
  assert_radius_below(acb_realref(MK1), 1e-12);
  assert_radius_below(acb_imagref(MK1), 1e-12);
  assert(arb_contains_zero(acb_imagref(MK1)));

  print_zeros(L, c);
  assert_zeros(L, c);
  printf("dedekind_class_number: %s OK\n", c.label);

  acb_clear(MK1);
  acb_clear(rhs);
  Lfunc_clear(L);
}

static const long POLY_3_1_23_1[] = {-1, -1, 0, 1}; // x^3 - x - 1
static const long BAD_3_1_23_1_23[] = {1, -1};
static const bad_factor BAD_3_1_23_1[] = {
  {23, BAD_3_1_23_1_23, 2},
};

static const long POLY_3_1_283_1[] = {-1, 4, 0, 1}; // x^3 + 4*x - 1
static const long BAD_3_1_283_1_283[] = {1, -1};
static const bad_factor BAD_3_1_283_1[] = {
  {283, BAD_3_1_283_1_283, 2},
};

static const long POLY_4_0_117_1[] = {1, 1, -1, -1, 1}; // x^4 - x^3 - x^2 + x + 1
static const long BAD_4_0_117_1_3[] = {1, 1};
static const long BAD_4_0_117_1_13[] = {1, 0, -1};
static const bad_factor BAD_4_0_117_1[] = {
  {3, BAD_4_0_117_1_3, 2},
  {13, BAD_4_0_117_1_13, 3},
};

static const long POLY_4_0_229_1[] = {1, -1, 0, 0, 1}; // x^4 - x + 1
static const long BAD_4_0_229_1_229[] = {1, 0, -1};
static const bad_factor BAD_4_0_229_1[] = {
  {229, BAD_4_0_229_1_229, 3},
};

static const long POLY_4_0_1521_1[] = {9, 3, 4, -1, 1}; // x^4 - x^3 + 4*x^2 + 3*x + 9
static const long BAD_4_0_1521_1_2[] = {1, 1, -1, -1};
static const long BAD_4_0_1521_1_3[] = {1, -1};
static const long BAD_4_0_1521_1_13[] = {1, -1};
static const bad_factor BAD_4_0_1521_1[] = {
  // 2 is an inessential prime for this non-monogenic defining polynomial:
  // the field is unramified at 2, but factorization mod 2 is not enough.
  {2, BAD_4_0_1521_1_2, 4},
  {3, BAD_4_0_1521_1_3, 2},
  {13, BAD_4_0_1521_1_13, 2},
};

static const long POLY_5_1_1609_1[] = {1, 1, -1, -1, 0, 1}; // x^5 - x^3 - x^2 + x + 1
static const long BAD_5_1_1609_1_1609[] = {1, -1, -1, 1};
static const bad_factor BAD_5_1_1609_1[] = {
  {1609, BAD_5_1_1609_1_1609, 4},
};

static const long POLY_6_0_14731_1[] = {1, 0, -1, 1, 0, -1, 1}; // x^6 - x^5 + x^3 - x^2 + 1
static const long BAD_6_0_14731_1_14731[] = {1, 0, 0, 0, -1};
static const bad_factor BAD_6_0_14731_1[] = {
  {14731, BAD_6_0_14731_1_14731, 5},
};

static const long POLY_7_1_184607_1[] = {1, 1, -1, 0, 1, -1, -1, 1};
static const long BAD_7_1_184607_1_184607[] = {1, -3, 2, 2, -3, 1};
static const bad_factor BAD_7_1_184607_1[] = {
  {184607, BAD_7_1_184607_1_184607, 6},
};

static const long POLY_8_0_1740113_1[] = {1, -1, 3, -4, 4, -4, 3, -2, 1};
static const long BAD_8_0_1740113_1_1740113[] = {1, 0, 0, -2, 0, 0, 1};
static const bad_factor BAD_8_0_1740113_1[] = {
  {1740113, BAD_8_0_1740113_1_1740113, 7},
};

static const dedekind_testcase DEDEKIND_TESTS[] = {
  // Cubic field 3.1.23.1.  The zero goldens come from PARI, with the
  // Riemann-zeta zero at 14.134725... removed:
  // lfunzeros(lfuncreate(x^3 - x - 1), 20).
  {
    "3.1.23.1",
    POLY_3_1_23_1, 4,
    BAD_3_1_23_1, 1,
    2, 23,
    {0.0, 1.0, 0.0, 0.0},
    1, 1,
    1,
    "0.28119957432296184651205076406787829979202322574406646267573",
    2, 23,
    0,
    "1", "0",
    {
      "5.11568332881511759855",
      "7.15926229054170384129",
      "8.88139657368917747496",
      "10.2820274008521320537",
      "11.4300363530481331189",
      "12.9344096677184118357",
      "14.6624809545746591259",
      "16.4982325133624013415",
      NULL,
    },
  },
  // Cubic field 3.1.283.1.  This exercises the class-number factor h_K = 2
  // while keeping the quotient degree 2.  Zero goldens come from PARI:
  // lfunzeros(lfuncreate(x^3 + 4*x - 1), 20).
  {
    "3.1.283.1",
    POLY_3_1_283_1, 4,
    BAD_3_1_283_1, 1,
    2, 283,
    {0.0, 1.0, 0.0, 0.0},
    1, 1,
    2,
    "1.401342327308812661850449007055876367036976238101096245011464371726157379243146620104485627027558603",
    2, 283,
    0,
    "1", "0",
    {
      "2.19353375194441133963",
      "3.43060415019741538455",
      "5.03445355437501277594",
      "6.31517765466829869223",
      "7.15341918498181596628",
      "7.93359343011215227188",
      "9.55723750264926372914",
      "10.1864841115568513881",
      NULL,
    },
  },
  // Quartic field 4.0.117.1.  This exercises a degree-3 quotient with two
  // ramified primes and w_K = 6.  Zero goldens come from PARI:
  // lfunzeros(lfuncreate(x^4 - x^3 - x^2 + x + 1), 20).
  {
    "4.0.117.1",
    POLY_4_0_117_1, 5,
    BAD_4_0_117_1, 2,
    3, 117,
    {0.0, 1.0, 1.0, 0.0},
    0, 2,
    1,
    "0.5435350724978695498926364006192337217181767537693045902578565450006843770349604030842035014653078917",
    6, 117,
    0,
    "1", "0",
    {
      "4.39582175376675801531",
      "5.72936210305730426417",
      "7.54499009544086055036",
      "8.03973715568146668171",
      "9.33426191622168877964",
      "10.4260852542660815246",
      "11.2492062077729352497",
      "11.9327135302546740372",
      NULL,
    },
  },
  // Quartic field 4.0.229.1, x^4 - x + 1.  Its quotient M_K is the same
  // degree-3 S4 Artin L-function tested in examples/artin.cpp as
  // 3.229.4t5.a.a, so the zero goldens below are copied from that vetted
  // Artin regression row.  They are independently reproduced by PARI:
  // lfunzeros(lfuncreate(x^4 - x + 1), 10).
  {
    "4.0.229.1",
    POLY_4_0_229_1, 5,
    BAD_4_0_229_1, 1,
    3, 229,
    {0.0, 1.0, 1.0, 0.0},
    0, 2,
    1,
    "0.3373778035715619014381982404044614505035293766460445399946133389725856960148538114465710607840662889",
    2, 229,
    0,
    "1", "0",
    {
      "3.19363569581669561149",
      "5.05331678840424750546",
      "6.46489851926996045565",
      "7.45209709798322283615",
      "8.62206109779680217600",
      "9.29879110967166199257",
      "10.0343434496103083426",
      "10.7698940399578528813",
      NULL,
    },
  },
  // Quartic field 4.0.1521.1.  This combines a degree-3 quotient with
  // h_K = 2 and two ramified primes.  Zero goldens come from PARI:
  // lfunzeros(lfuncreate(x^4 - x^3 + 4*x^2 + 3*x + 9), 20).
  {
    "4.0.1521.1",
    POLY_4_0_1521_1, 5,
    BAD_4_0_1521_1, 3,
    3, 1521,
    {0.0, 1.0, 1.0, 0.0},
    0, 2,
    2,
    "2.389526434574218608223861657038181047072324150306010858541360598922648191661925060537743228578628751",
    6, 1521,
    0,
    "1", "0",
    {
      "2.15367543022529411299",
      "3.11934147900860341390",
      "3.95537331939419894859",
      "5.85202180427453258469",
      "6.79945605063326105935",
      "7.23159073941876201502",
      "8.03973715568146668171",
      "8.62542663503259159734",
      NULL,
    },
  },
  // Quintic field 5.1.1609.1.  This is the first quotient-degree-4 row:
  // signature (1,2), h_K = 1, w_K = 2, and one ramified prime.  Zero goldens
  // come from PARI:
  // lfunzeros(lfuncreate(x^5 - x^3 - x^2 + x + 1), 20).
  {
    "5.1.1609.1",
    POLY_5_1_1609_1, 6,
    BAD_5_1_1609_1, 1,
    4, 1609,
    {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
    1, 2,
    1,
    "0.2683555508375814276911063300754518114345280574289672119596266683595315375418520115002384600100240896",
    2, 1609,
    0,
    "1", "0",
    {
      "3.50464340448417603999",
      "4.40290813303232740396",
      "5.85582736823442197038",
      "6.37980089698034500761",
      "7.13618003079659589336",
      "8.24029657577021448469",
      "8.72066903279375410192",
      "9.33184883095975093865",
      NULL,
    },
  },
  // Sextic field 6.0.14731.1.  This gives a quotient-degree-5 row:
  // signature (0,3), h_K = 1, w_K = 2, and one ramified prime.  Zero goldens
  // come from PARI:
  // lfunzeros(lfuncreate(x^6 - x^5 + x^3 - x^2 + 1), 15).
  {
    "6.0.14731.1",
    POLY_6_0_14731_1, 7,
    BAD_6_0_14731_1, 1,
    5, 14731,
    {0.0, 0.0, 1.0, 1.0, 1.0, 0.0},
    0, 3,
    1,
    "0.2777544599776168607506758037114708999530031988686870787222596040344583176156311231176770855555381710",
    2, 14731,
    0,
    "1", "0",
    {
      "2.76421217155354135232",
      "3.85705987401773225118",
      "4.96457136377154903045",
      "5.52884541100029020892",
      "6.16311996426779467920",
      "6.81431539802753313405",
      "7.40240817630714278728",
      "7.97221102266334064496",
      NULL,
    },
  },
  // Septic field 7.1.184607.1.  This gives a quotient-degree-6 row:
  // signature (1,3), h_K = 1, w_K = 2, and one ramified prime.  Zero goldens
  // come from PARI:
  // lfunzeros(lfuncreate(x^7 - x^6 - x^5 + x^4 - x^2 + x + 1), 12).
  {
    "7.1.184607.1",
    POLY_7_1_184607_1, 8,
    BAD_7_1_184607_1, 1,
    6, 184607,
    {0.0, 0.0, 0.0, 1.0, 1.0, 1.0},
    1, 3,
    1,
    "0.3804471063197960076887843706074230749712618217711544276546078797203530270399830051145265834434539766",
    2, 184607,
    0,
    "1", "0",
    {
      "2.90019611804893278263",
      "3.22406920147009107702",
      "4.42828461462662743355",
      "4.82495879790651215688",
      "5.51960702128344437712",
      "5.86860611718421182072",
      "6.62669377214100435226",
      "7.11526553535783770214",
      NULL,
    },
  },
  // Degree-8 field 8.0.1740113.1.  This gives a quotient-degree-7 stress
  // case.  Zero goldens come from PARI with an 8 GB stack:
  // allocatemem(8000000000);
  // lfunzeros(lfuncreate(x^8 - 2*x^7 + 3*x^6 - 4*x^5 + 4*x^4 - 4*x^3
  //                      + 3*x^2 - x + 1), 12).
  {
    "8.0.1740113.1",
    POLY_8_0_1740113_1, 9,
    BAD_8_0_1740113_1, 1,
    7, 1740113,
    {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0},
    0, 4,
    1,
    "0.3887666677436626597100277928686136571894987881092570754354061379605713862101060251365705063515929098",
    2, 1740113,
    0,
    "1", "0",
    {
      "2.36339922245404673570",
      "3.21706680342524374400",
      "3.80064136777277898728",
      "4.53189984191840889327",
      "4.91089046416537679892",
      "5.34075844793792231447",
      "5.83809561853187865141",
      "6.42140405113969243841",
      NULL,
    },
  },
};

int main()
{
  assert_lpoly_sanity_checks();

  const size_t ncases = sizeof(DEDEKIND_TESTS) / sizeof(DEDEKIND_TESTS[0]);
  for (size_t i = 0; i < ncases; ++i)
    run_case(DEDEKIND_TESTS[i]);

  flint_cleanup();
  return 0;
}
