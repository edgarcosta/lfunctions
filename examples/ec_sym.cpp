// Copyright Edgar Costa 2026
// See LICENSE file for license details.
/*
 * Worked symmetric-power L-functions of the elliptic curve 11.a1 (bead lfunctions-cu2).
 *
 * Builds the Sym^2 (degree 3) and Sym^3 (degree 4) L-functions of
 * 11.a1 = [0,-1,1,-10,-20] and asserts their certified rank, root number, first
 * zero, and an off-central special value.
 *
 *  - Good-prime Euler factors are formed in C from the curve's a_p by
 *    sym_power_lpoly (bead lfunctions-4q6).  The a_p are hardcoded (the library
 *    has no elliptic-curve point counting); the table was generated once with
 *    smalljac's lpdata and matches the LMFDB newform 11.2.a.a.
 *  - 11.a1 has split multiplicative reduction at 11 (a_11 = +1), so the only bad
 *    factor is 1 - a_11^k T = 1 - T = [1,-1] at p = 11 for every k.  This
 *    deliberately avoids the additive-prime trap: additive/ramified bad factors
 *    are not a function of a_p and so are out of scope for the helper.
 *  - self_dual = YES (a symmetric power of the self-dual L-function of an
 *    elliptic curve is self-dual), declared through the public Lfunc_init_advanced
 *    rather than poking the opaque handle, mirroring test/highdeg_check.cpp.
 *  - Uses the default Buthe verifier, which certifies RH-completeness for
 *    degrees 2-9; any ERR_RH_ERROR would be a non-fatal warning, so it is
 *    collected and reported, never asserted against.
 *
 * Goldens: Sym^2 rank / eps / first zero from the certified high-degree suite
 * (test/highdeg/objects.yaml); Sym^3 from Pari lfunsympow(E,3) via lfunorderzero
 * / lfunrootres / lfunzeros; both L(2.5) values from Pari lfun, in the algebraic
 * (arithmetic) normalisation that Lfunc_special_value and Pari lfun share.  The
 * central leading-Taylor coefficient is deliberately NOT asserted: its sympow
 * normalisation differs between Pari and the library (objects.yaml omits taylor
 * for sympow), so an unambiguous off-central value is used instead.
 *
 * Expected output (assertions, not eyeballed -- a failure aborts nonzero):
 *   Sym^2 of 11.a1 (degree 3, conductor 121):  rank 0, eps (1,0),
 *     first zero 3.8992814947713447822425023053642116570, L(2.5) 1.1197003417627197831
 *   Sym^3 of 11.a1 (degree 4, conductor 1331): rank 0, eps (1,0),
 *     first zero 2.3200209723548958126909115564099704820, L(2.5) 1.1763280758703274817
 */
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "sym_power.h"
#include <cassert>

// a_p of 11.a1 = [0,-1,1,-10,-20] at every good prime p <= 5065 (= the Sym^3
// nmax; Sym^2 needs only p <= 692).  Generated once with
// `lpdata ap "[0,-1,1,-10,-20]" 5065 1` (a_p = -t); the bad prime 11 is
// intentionally absent.  Matches LMFDB newform 11.2.a.a.
static const std::map<uint64_t, int64_t> ap_11a1 = {
  {2, -2}, {3, -1}, {5, 1}, {7, -2}, {13, 4}, {17, -2}, {19, 0}, {23, -1},
  {29, 0}, {31, 7}, {37, 3}, {41, -8}, {43, -6}, {47, 8}, {53, -6}, {59, 5},
  {61, 12}, {67, -7}, {71, -3}, {73, 4}, {79, -10}, {83, -6}, {89, 15}, {97, -7},
  {101, 2}, {103, -16}, {107, 18}, {109, 10}, {113, 9}, {127, 8}, {131, -18}, {137, -7},
  {139, 10}, {149, -10}, {151, 2}, {157, -7}, {163, 4}, {167, -12}, {173, -6}, {179, -15},
  {181, 7}, {191, 17}, {193, 4}, {197, -2}, {199, 0}, {211, 12}, {223, 19}, {227, 18},
  {229, 15}, {233, 24}, {239, -30}, {241, -8}, {251, -23}, {257, -2}, {263, 14}, {269, 10},
  {271, -28}, {277, -2}, {281, -18}, {283, 4}, {293, 24}, {307, 8}, {311, 12}, {313, -1},
  {317, 13}, {331, 7}, {337, -22}, {347, 28}, {349, 30}, {353, -21}, {359, -20}, {367, -17},
  {373, -26}, {379, -5}, {383, -1}, {389, -15}, {397, -2}, {401, 2}, {409, -30}, {419, 20},
  {421, 22}, {431, -18}, {433, -11}, {439, 40}, {443, -11}, {449, 35}, {457, -12}, {461, 12},
  {463, -11}, {467, -27}, {479, 20}, {487, 23}, {491, -8}, {499, 20}, {503, -26}, {509, 15},
  {521, -3}, {523, -16}, {541, -8}, {547, 8}, {557, -2}, {563, 4}, {569, 0}, {571, -28},
  {577, 33}, {587, 28}, {593, 44}, {599, 40}, {601, 2}, {607, -22}, {613, -16}, {617, 18},
  {619, -25}, {631, 7}, {641, -33}, {643, 29}, {647, -7}, {653, -41}, {659, 10}, {661, 37},
  {673, 14}, {677, -42}, {683, -16}, {691, 17}, {701, 2}, {709, -25}, {719, 15}, {727, 3},
  {733, -36}, {739, 50}, {743, 4}, {751, -23}, {757, -22}, {761, 12}, {769, 20}, {773, -6},
  {787, -32}, {797, 53}, {809, 0}, {811, -38}, {821, 22}, {823, 39}, {827, -52}, {829, 25},
  {839, -5}, {853, 14}, {857, 8}, {859, -15}, {863, 24}, {877, -12}, {881, -43}, {883, 4},
  {887, -22}, {907, -12}, {911, 12}, {919, 10}, {929, -30}, {937, 8}, {941, 42}, {947, -27},
  {953, 34}, {967, -32}, {971, 47}, {977, -27}, {983, 39}, {991, -8}, {997, 38}, {1009, -10},
  {1013, 39}, {1019, -10}, {1021, 22}, {1031, 32}, {1033, -16}, {1039, 5}, {1049, -55}, {1051, 2},
  {1061, -13}, {1063, 44}, {1069, -20}, {1087, 8}, {1091, -58}, {1093, -51}, {1097, -42}, {1103, -51},
  {1109, -30}, {1117, 48}, {1123, 24}, {1129, 50}, {1151, 2}, {1153, -31}, {1163, 34}, {1171, -3},
  {1181, -18}, {1187, -12}, {1193, -21}, {1201, 2}, {1213, -41}, {1217, -42}, {1223, 14}, {1229, 60},
  {1231, -18}, {1237, 18}, {1249, 40}, {1259, -25}, {1277, -47}, {1279, -15}, {1283, -36}, {1289, 0},
  {1291, -8}, {1297, 48}, {1301, 27}, {1303, 39}, {1307, 28}, {1319, -30}, {1321, 47}, {1327, 68},
  {1361, 12}, {1367, -72}, {1373, 39}, {1381, -68}, {1399, 60}, {1409, -15}, {1423, 29}, {1427, -12},
  {1429, -70}, {1433, 54}, {1439, 0}, {1447, 28}, {1451, 52}, {1453, -71}, {1459, -20}, {1471, 22},
  {1481, 32}, {1483, 49}, {1487, 58}, {1489, -15}, {1493, -36}, {1499, 55}, {1511, 37}, {1523, -41},
  {1531, 32}, {1543, -36}, {1549, -15}, {1553, -56}, {1559, -60}, {1567, -52}, {1571, -28}, {1579, -30},
  {1583, 34}, {1597, -32}, {1601, 2}, {1607, 33}, {1609, -10}, {1613, -6}, {1619, -20}, {1621, 22},
  {1627, 78}, {1637, 33}, {1657, -2}, {1663, 4}, {1667, 48}, {1669, 50}, {1693, -6}, {1697, -42},
  {1699, 40}, {1709, -45}, {1721, -3}, {1723, -46}, {1733, -6}, {1741, 17}, {1747, -57}, {1753, 34},
  {1759, -40}, {1777, 8}, {1783, 59}, {1787, -57}, {1789, 10}, {1801, 52}, {1811, 12}, {1823, -56},
  {1831, -43}, {1847, -52}, {1861, 62}, {1867, 28}, {1871, -3}, {1873, -6}, {1877, 18}, {1879, -35},
  {1889, 70}, {1901, 77}, {1907, -52}, {1913, -36}, {1931, -18}, {1933, 54}, {1949, -40}, {1951, -23},
  {1973, 79}, {1979, 30}, {1987, -22}, {1993, -66}, {1997, -72}, {1999, -20}, {2003, 4}, {2011, -13},
  {2017, -17}, {2027, 63}, {2029, 45}, {2039, 60}, {2053, 84}, {2063, 24}, {2069, 70}, {2081, -18},
  {2083, 89}, {2087, 48}, {2089, -10}, {2099, 35}, {2111, -38}, {2113, -86}, {2129, 20}, {2131, -68},
  {2137, 73}, {2141, -58}, {2143, -91}, {2153, -26}, {2161, -13}, {2179, -45}, {2203, -1}, {2207, 48},
  {2213, 4}, {2221, 22}, {2237, 78}, {2239, -70}, {2243, -56}, {2251, -48}, {2267, 93}, {2269, 25},
  {2273, 4}, {2281, 7}, {2287, 38}, {2293, 29}, {2297, -57}, {2309, 60}, {2311, -13}, {2333, 59},
  {2339, 10}, {2341, 67}, {2347, -37}, {2351, -48}, {2357, -57}, {2371, -28}, {2377, 3}, {2381, -18},
  {2383, -36}, {2389, -50}, {2393, 54}, {2399, -75}, {2411, 62}, {2417, -22}, {2423, -31}, {2437, -82},
  {2441, 42}, {2447, 3}, {2459, -50}, {2467, 3}, {2473, -11}, {2477, 48}, {2503, 14}, {2521, 72},
  {2531, 57}, {2539, 0}, {2543, 34}, {2549, -20}, {2551, -98}, {2557, 13}, {2579, 20}, {2591, -58},
  {2593, 14}, {2609, -30}, {2617, 18}, {2621, 22}, {2633, 39}, {2647, 38}, {2657, 38}, {2659, 40},
  {2663, 39}, {2671, 72}, {2677, -7}, {2683, -16}, {2687, 23}, {2689, 5}, {2693, -41}, {2699, -55},
  {2707, -17}, {2711, 87}, {2713, -56}, {2719, -70}, {2729, 30}, {2731, -68}, {2741, -58}, {2749, 50},
  {2753, 49}, {2767, 48}, {2777, -42}, {2789, -20}, {2791, 42}, {2797, -42}, {2801, 52}, {2803, 44},
  {2819, -25}, {2833, -6}, {2837, -62}, {2843, 4}, {2851, 2}, {2857, -82}, {2861, -63}, {2879, -40},
  {2887, -57}, {2897, 38}, {2903, 54}, {2909, -25}, {2917, 88}, {2927, -72}, {2939, -50}, {2953, -86},
  {2957, 3}, {2963, -81}, {2969, 70}, {2971, -53}, {2999, -80}, {3001, 27}, {3011, 62}, {3019, 85},
  {3023, 39}, {3037, 13}, {3041, 42}, {3049, -40}, {3061, 37}, {3067, 13}, {3079, -20}, {3083, 29},
  {3089, 25}, {3109, 80}, {3119, -90}, {3121, 22}, {3137, 8}, {3163, -26}, {3167, 18}, {3169, 45},
  {3181, 32}, {3187, -2}, {3191, -8}, {3203, -6}, {3209, -10}, {3217, 23}, {3221, -103}, {3229, 70},
  {3251, -48}, {3253, 74}, {3257, 58}, {3259, -60}, {3271, -3}, {3299, 100}, {3301, -73}, {3307, 98},
  {3313, 4}, {3319, 0}, {3323, -76}, {3329, -100}, {3331, -43}, {3343, 44}, {3347, -17}, {3359, -45},
  {3361, -88}, {3371, -103}, {3373, 4}, {3389, 15}, {3391, 92}, {3407, 18}, {3413, 9}, {3433, -66},
  {3449, 40}, {3457, -57}, {3461, -38}, {3463, -96}, {3467, 38}, {3469, -85}, {3491, 17}, {3499, 100},
  {3511, 12}, {3517, -22}, {3527, 18}, {3529, -35}, {3533, 24}, {3539, 20}, {3541, 42}, {3547, 53},
  {3557, -27}, {3559, 0}, {3571, -28}, {3581, 32}, {3583, -96}, {3593, -26}, {3607, 58}, {3613, -26},
  {3617, 63}, {3623, 4}, {3631, 32}, {3637, -72}, {3643, 34}, {3659, 30}, {3671, -78}, {3673, -76},
  {3677, -62}, {3691, 92}, {3697, 23}, {3701, 102}, {3709, 20}, {3719, -55}, {3727, 23}, {3733, 114},
  {3739, -60}, {3761, -88}, {3767, -27}, {3769, 40}, {3779, -30}, {3793, 34}, {3797, -82}, {3803, 74},
  {3821, -3}, {3823, 84}, {3833, 19}, {3847, -42}, {3851, -73}, {3853, 74}, {3863, 54}, {3877, 58},
  {3881, 7}, {3889, -70}, {3907, -22}, {3911, 12}, {3917, -57}, {3919, 0}, {3923, 54}, {3929, -60},
  {3931, 107}, {3943, -96}, {3947, -107}, {3967, -92}, {3989, -90}, {4001, 102}, {4003, -46}, {4007, -32},
  {4013, 54}, {4019, 15}, {4021, 22}, {4027, 28}, {4049, 25}, {4051, -123}, {4057, 103}, {4073, -31},
  {4079, -60}, {4091, 92}, {4093, 94}, {4099, 20}, {4111, 62}, {4127, 48}, {4129, -25}, {4133, 4},
  {4139, 20}, {4153, -36}, {4157, -42}, {4159, 85}, {4177, 88}, {4201, -98}, {4211, -63}, {4217, 33},
  {4219, 10}, {4229, -55}, {4231, -68}, {4241, 92}, {4243, 64}, {4253, -16}, {4259, 70}, {4261, -113},
  {4271, -53}, {4273, -41}, {4283, 39}, {4289, -60}, {4297, -62}, {4327, -107}, {4337, -2}, {4339, 55},
  {4349, -30}, {4357, -117}, {4363, -6}, {4373, 84}, {4391, 42}, {4397, 108}, {4409, -30}, {4421, 72},
  {4423, 49}, {4441, 42}, {4447, 83}, {4451, 102}, {4457, 78}, {4463, -66}, {4481, 82}, {4483, -106},
  {4493, 29}, {4507, -82}, {4513, 59}, {4517, 8}, {4519, 5}, {4523, -36}, {4547, -12}, {4549, -10},
  {4561, 112}, {4567, -112}, {4583, -36}, {4591, -8}, {4597, -52}, {4603, 89}, {4621, 122}, {4637, 18},
  {4639, 20}, {4643, -81}, {4649, -80}, {4651, -73}, {4657, 73}, {4663, 64}, {4673, 114}, {4679, 25},
  {4691, 17}, {4703, -36}, {4721, 22}, {4723, -96}, {4729, 30}, {4733, 129}, {4751, -48}, {4759, 30},
  {4783, -101}, {4787, 8}, {4789, 110}, {4793, -36}, {4799, 105}, {4801, 77}, {4813, -76}, {4817, -132},
  {4831, -68}, {4861, -88}, {4871, 72}, {4877, 93}, {4889, 5}, {4903, -126}, {4909, 115}, {4919, 30},
  {4931, 32}, {4933, 94}, {4937, -42}, {4943, -76}, {4951, -48}, {4957, 98}, {4967, -52}, {4969, 100},
  {4973, -1}, {4987, 28}, {4993, 44}, {4999, 40}, {5003, 119}, {5009, -45}, {5011, -38}, {5021, -3},
  {5023, 54}, {5039, 15}, {5051, -48}, {5059, 10},
};

// Which symmetric power the lpoly callback should form.
struct sym_ctx { int k; };

// Euler factor of Sym^k(11.a1) at p.  At the split-multiplicative bad prime 11
// it is 1 - a_11^k T = 1 - T (a_11 = +1) for every k; at good primes it is
// sym_power_lpoly(a_p, p, k) (exact fmpz) converted to acb at working precision.
static void sym_lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)),
                               int64_t prec, void *param)
{
  const int k = static_cast<sym_ctx *>(param)->k;
  acb_poly_zero(poly);
  if (p == 11) {
    acb_poly_set_coeff_si(poly, 0, 1);
    acb_poly_set_coeff_si(poly, 1, -1);
    return;
  }
  const auto it = ap_11a1.find(p);
  if (it == ap_11a1.end()) {
    // A good prime <= nmax missing from the table would silently corrupt the
    // result: returning the zero poly makes the library OR in ERR_INSUFF_EULER,
    // which is only a warning, not fatal.  Fail loudly even when assertions are
    // compiled out.  The table is sized to the current Sym^2/Sym^3 nmax
    // (692/5065); a larger conductor or precision could outgrow it.
    fprintf(stderr, "ec_sym: a_p table missing good prime %lu (<= nmax)\n", (unsigned long) p);
    std::abort();
  }
  fmpz_poly_t f;
  fmpz_poly_init(f);
  sym_power_lpoly(f, it->second, p, k);
  acb_poly_set_fmpz_poly(poly, f, prec);
  fmpz_poly_clear(f);
}

// Build Sym^k of 11.a1 and assert its certified rank / |eps| = 1 / eps / first
// zero, plus the off-central special value L(2.5) against its Pari golden.
// Returns accumulated warnings (notably the tolerated ERR_RH_ERROR for degree
// >= 3); a fatal error or a failed assertion aborts.
static Lerror_t check_sym_power(int k, uint64_t conductor, double normalisation,
                                double *mus, int expected_rank, double eps_re, double eps_im,
                                const char *first_zero, const char *l_at_2p5)
{
  const slong prec = 400;
  Lerror_t ecode = 0;
  char cache_dir[] = ".";

  Lparams_t Lp;
  Lp.degree = (uint64_t) (k + 1);
  Lp.conductor = conductor;
  Lp.normalisation = normalisation;
  Lp.mus = mus;                 // Lfunc_init_advanced copies this
  Lp.target_prec = DEFAULT_TARGET_PREC;
  Lp.wprec = 0;
  Lp.gprec = 0;
  Lp.self_dual = YES;           // Sym^k of a self-dual EC L-function is self-dual
  Lp.rank = DK;
  Lp.cache_dir = cache_dir;

  Lfunc_t L = Lfunc_init_advanced(&Lp, &ecode);
  assert(!fatal_error(ecode));

  sym_ctx ctx;
  ctx.k = k;
  ecode |= Lfunc_use_all_lpolys(L, sym_lpoly_callback, &ctx);
  assert(!fatal_error(ecode));

  ecode |= Lfunc_compute(L);    // any ERR_RH_ERROR is a non-fatal warning, collected not asserted
  assert(!fatal_error(ecode));

  printf("Sym^%d of 11.a1 (degree %d, conductor %lu):\n", k, k + 1, (unsigned long) conductor);

  // (1) rank
  printf("  rank       = %ld (want %d)\n", (long) Lfunc_rank(L), expected_rank);
  assert(Lfunc_rank(L) == expected_rank);

  // (2) |eps| = 1, and eps == (eps_re, eps_im); each reference widened by 1e-9.
  arb_t mag, one, err;
  arb_init(mag); arb_init(one); arb_init(err);
  acb_abs(mag, Lfunc_sign(L), prec);
  arb_set_ui(one, 1);
  arb_set_d(err, 1e-9); arb_add_error(one, err);
  assert(arb_overlaps(mag, one));
  acb_t epsref;
  acb_init(epsref);
  acb_set_d_d(epsref, eps_re, eps_im);
  arb_set_d(err, 1e-9); acb_add_error_arb(epsref, err);
  printf("  eps        = "); acb_printd(Lfunc_sign(L), 20); printf("  (want (%g, %g))\n", eps_re, eps_im);
  assert(acb_overlaps(Lfunc_sign(L), epsref));
  arb_clear(mag); arb_clear(one); arb_clear(err); acb_clear(epsref);

  // (3) first zero, against [first_zero +- 1e-15].
  arb_t zref, zerr;
  arb_init(zref); arb_init(zerr);
  arb_set_str(zref, first_zero, prec);
  arb_set_d(zerr, 1e-15); arb_add_error(zref, zerr);
  printf("  first zero = "); arb_printd(Lfunc_zeros(L, 0) + 0, 30);
  printf("\n               (want %s)\n", first_zero);
  assert(arb_overlaps(Lfunc_zeros(L, 0) + 0, zref));
  arb_clear(zref); arb_clear(zerr);

  // (4) off-central special value L(2.5) (real), against the Pari golden, with
  //     absolute tolerance 1e-12.  Lfunc_special_value takes its argument in the
  //     algebraic (arithmetic) normalisation -- the same one Pari lfun uses for
  //     lfunsympow -- so 2.5 maps directly (Sym^k is centred at (k+1)/2; 2.5 is
  //     off-centre and real, where the value is unambiguous).
  acb_t lv, lref;
  acb_init(lv); acb_init(lref);
  ecode |= Lfunc_special_value(lv, L, 2.5, 0.0);
  assert(!fatal_error(ecode));
  arb_set_str(acb_realref(lref), l_at_2p5, prec);
  arb_set_str(acb_imagref(lref), "0", prec);
  arb_t lerr; arb_init(lerr); arb_set_d(lerr, 1e-12); acb_add_error_arb(lref, lerr); arb_clear(lerr);
  printf("  L(2.5)     = "); acb_printd(lv, 20); printf("  (want %s)\n", l_at_2p5);
  assert(acb_overlaps(lv, lref));
  acb_clear(lv); acb_clear(lref);

  Lfunc_clear(L);
  return ecode;
}

int main(void)
{
  Lerror_t ecode = 0;

  // Sym^2: degree 3, conductor 121, norm 1.0, mus [0,1,0].  rank/eps/first-zero
  // golden from the high-degree suite (test/highdeg/objects.yaml); L(2.5) from Pari.
  double mus2[3] = {0, 1, 0};
  ecode |= check_sym_power(2, 121, 1.0, mus2, 0, 1.0, 0.0,
                           "3.8992814947713447822425023053642116570",
                           "1.119700341762719783127772146865354330659");

  // Sym^3: degree 4, conductor 1331, norm 1.5, mus [0,-1,1,0].  Golden from Pari
  // lfunsympow(E,3): lfunorderzero / lfunrootres / lfunzeros, and lfun for L(2.5).
  double mus3[4] = {0, -1, 1, 0};
  ecode |= check_sym_power(3, 1331, 1.5, mus3, 0, 1.0, 0.0,
                           "2.3200209723548958126909115564099704819803465711244",
                           "1.176328075870327481688239970219759609463");

  fprint_errors(stderr, ecode);   // report collected warnings (tolerated ERR_RH_ERROR)
  printf("ec_sym: all assertions passed\n");
  return 0;
}
