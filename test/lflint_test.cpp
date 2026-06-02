// Unit tests for lfun::ZZ and lfun::ZZX (lflint.h).
// Plain <cassert>; exits 0 on success.
#include <cassert>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include "lflint.h"

using lfun::ZZ;
using lfun::ZZX;

static void test_default_ctor_is_zero() {
  ZZ a;
  assert(a.is_zero());
}

static void test_construction_and_raw() {
  // ctor from slong / ulong
  ZZ a(slong{42});
  ZZ b(ulong{7});
  assert(!a.is_zero() && !b.is_zero());
  // raw accessors
  assert(fmpz_get_si(a._fmpz()) == 42);
  assert(fmpz_get_ui(b._fmpz()) == 7);

  // copy ctor
  ZZ c(a);
  assert(fmpz_equal(c._fmpz(), a._fmpz()));

  // move ctor leaves source destructible (== 0 by design)
  ZZ d(std::move(c));
  assert(fmpz_get_si(d._fmpz()) == 42);
  assert(c.is_zero());

  // copy-assign
  ZZ e;
  e = a;
  assert(fmpz_equal(e._fmpz(), a._fmpz()));

  // move-assign
  ZZ f;
  f = std::move(e);
  assert(fmpz_get_si(f._fmpz()) == 42);
  assert(e.is_zero());

  // assign from slong / ulong
  ZZ g;
  g = slong{-13};
  assert(fmpz_get_si(g._fmpz()) == -13);
  g = ulong{9};
  assert(fmpz_get_ui(g._fmpz()) == 9);
}

static void test_set_zero_one() {
  ZZ a(slong{42});
  a.set_zero();
  assert(a.is_zero());
  a.set_one();
  assert(a.is_one());
  assert(!a.is_zero());
}

static void test_arithmetic() {
  ZZ two(slong{2});
  ZZ three(slong{3});

  // unary -
  ZZ neg = -two;
  assert(fmpz_get_si(neg._fmpz()) == -2);

  // binary +, -, *
  ZZ sum  = two + three;
  ZZ diff = three - two;
  ZZ prod = two * three;
  assert(fmpz_get_si(sum._fmpz())  == 5);
  assert(fmpz_get_si(diff._fmpz()) == 1);
  assert(fmpz_get_si(prod._fmpz()) == 6);

  // slong * ZZ and ZZ * slong (used in sympow_ECQ: `9*a[8]`, `a[2]*p[1]`)
  ZZ left  = slong{9} * three;
  ZZ right = three * slong{9};
  assert(fmpz_get_si(left._fmpz())  == 27);
  assert(fmpz_get_si(right._fmpz()) == 27);

  // compound
  ZZ x(slong{10});
  x += two;
  assert(fmpz_get_si(x._fmpz()) == 12);
  x *= three;
  assert(fmpz_get_si(x._fmpz()) == 36);
  x *= slong{2};
  assert(fmpz_get_si(x._fmpz()) == 72);

  // ZZ *= ulong (fmpz_mul_ui) — multiply by an unsigned value without
  // narrowing it through slong (pattern `qn *= q`, q a ulong)
  ZZ y(slong{6});
  y *= ulong{7};
  assert(fmpz_get_si(y._fmpz()) == 42);
}

static void test_modulo() {
  ZZ neg(slong{-5});
  ZZ three(slong{3});

  // ZZ % ZZ — fmpz_mod, non-negative
  ZZ r1 = neg % three;
  assert(fmpz_get_si(r1._fmpz()) == 1);   // (-5) mod 3 == 1 (NOT -2)

  // ZZ % ulong — fmpz_mod_ui, non-negative
  ZZ r2 = neg % ulong{3};
  assert(fmpz_get_si(r2._fmpz()) == 1);

  // positive sanity
  ZZ ten(slong{10});
  assert(fmpz_get_si((ten % three)._fmpz()) == 1);
  assert(fmpz_get_si((ten % ulong{3})._fmpz()) == 1);
}

static void test_comparison_and_map() {
  ZZ a(slong{5});
  ZZ b(slong{5});
  ZZ c(slong{7});

  assert(a == b);
  assert(!(a == c));
  assert(a != c);
  assert(!(a != b));
  assert(a == slong{5});
  assert(a != slong{6});

  // ordering (for std::map keys)
  assert(a < c);
  assert(!(c < a));
  assert(!(a < b));    // equal — strict less must be false

  // actual std::map usage
  std::map<ZZ, int> m;
  m[ZZ(slong{2})] = 20;
  m[ZZ(slong{1})] = 10;
  m[ZZ(slong{3})] = 30;
  auto it = m.begin();
  assert(it->second == 10); ++it;
  assert(it->second == 20); ++it;
  assert(it->second == 30);
}

static void test_divisibility() {
  ZZ minus_twelve(slong{-12});
  assert( minus_twelve.divisible(slong{4}));
  assert(!minus_twelve.divisible(slong{5}));
  ZZ q = minus_twelve.divexact(slong{4});
  assert(fmpz_get_si(q._fmpz()) == -3);

  // realistic usage — strip all factors of 4: `while (d.divisible(4)) d = d.divexact(4);`
  ZZ d(slong{-48});
  while (d.divisible(slong{4})) d = d.divexact(slong{4});
  assert(fmpz_get_si(d._fmpz()) == -3);
}

static void test_pow() {
  ZZ two(slong{2});
  ZZ p10 = lfun::pow(two, ulong{10});
  assert(fmpz_get_si(p10._fmpz()) == 1024);

  ZZ minus_three(slong{-3});
  ZZ p3 = lfun::pow(minus_three, ulong{3});
  assert(fmpz_get_si(p3._fmpz()) == -27);

  ZZ p0 = lfun::pow(minus_three, ulong{0});
  assert(p0.is_one());
}

static void test_to_conversion() {
  ZZ x(slong{-42});
  ZZ y(slong{1024});
  ZZ u(ulong{7});

  assert(x.to<slong>()  == -42);
  assert(y.to<slong>()  == 1024);
  assert(u.to<ulong>()  == ulong{7});
  assert(y.to<double>() == 1024.0);
  // realistic usage — reduce mod n, then narrow: `(coeff % n).to<mp_limb_t>()`
  // mp_limb_t == ulong; same specialisation
  assert((y % ulong{17}).to<ulong>() == ulong{1024 % 17});
}

static void test_ostream() {
  using std::ostringstream;
  using std::string;

  ostringstream os1; os1 << ZZ(slong{42});
  assert(os1.str() == string("42"));

  ostringstream os2; os2 << ZZ(slong{-13});
  assert(os2.str() == string("-13"));

  // multi-limb (10^25 doesn't fit in 64 bits)
  ZZ big;
  fmpz_set_str(big._fmpz(), "12345678901234567890123", 10);
  ostringstream os3; os3 << big;
  assert(os3.str() == string("12345678901234567890123"));
}

static void test_zzx_basic() {
  ZZX p;
  assert(p.degree() == -1);          // empty poly has degree -1 per FLINT

  ZZX q(slong{8});                   // capacity hint; still degree -1
  assert(q.degree() == -1);

  // set/get round-trip including negative coeff
  p.set_coeff(slong{3}, slong{-5});
  p.set_coeff(slong{1}, ZZ(slong{7}));
  assert(p.degree() == 3);
  assert(p.get_coeff(slong{3}).to<slong>() == -5);
  assert(p.get_coeff(slong{1}).to<slong>() ==  7);
  assert(p.get_coeff(slong{2}).is_zero());
  assert(p.get_coeff(slong{0}).is_zero());

  // fit_length does not raise degree if no coefficients beyond it are set
  ZZX r;
  r.set_coeff(slong{0}, slong{1});
  r.fit_length(slong{16});
  assert(r.degree() == 0);

  // operator= from int (constant) — used in examples_tools.h: `f = 0;`
  ZZX s;
  s.set_coeff(slong{2}, slong{99});
  assert(s.degree() == 2);
  s = slong{0};
  assert(s.degree() == -1);

  // copy ctor
  ZZX t(p);
  assert(t.get_coeff(slong{3}).to<slong>() == -5);
  assert(t.degree() == p.degree());

  // ZZX move ctor leaves source in degree -1 (valid, destructible)
  ZZX donor(p);                       // copy of p; degree 3
  ZZX recv(std::move(donor));
  assert(recv.degree() == 3);
  assert(donor.degree() == -1);

  // ZZX move-assign: swap semantics means source gets the target's old state
  ZZX target;                         // initially degree -1 (empty)
  target = std::move(recv);           // swap: target ← recv's state (degree 3), recv ← target's old state (degree -1)
  assert(target.degree() == 3);
  assert(recv.degree() == -1);
}

static void test_zzx_arithmetic() {
  // p = 1 + 2x + 3x^2
  ZZX p;
  p.set_coeff(slong{0}, slong{1});
  p.set_coeff(slong{1}, slong{2});
  p.set_coeff(slong{2}, slong{3});

  // q = 4 + 5x
  ZZX q;
  q.set_coeff(slong{0}, slong{4});
  q.set_coeff(slong{1}, slong{5});

  // p += q -> 5 + 7x + 3x^2
  ZZX r(p);
  r += q;
  assert(r.degree() == 2);
  assert(r.get_coeff(slong{0}).to<slong>() == 5);
  assert(r.get_coeff(slong{1}).to<slong>() == 7);
  assert(r.get_coeff(slong{2}).to<slong>() == 3);

  // (1 + T) * (1 - T) == 1 - T^2 — covers degree-shrink via cancellation
  ZZX one_plus_T;
  one_plus_T.set_coeff(slong{0}, slong{1});
  one_plus_T.set_coeff(slong{1}, slong{1});
  ZZX one_minus_T;
  one_minus_T.set_coeff(slong{0}, slong{ 1});
  one_minus_T.set_coeff(slong{1}, slong{-1});
  one_plus_T *= one_minus_T;
  assert(one_plus_T.degree() == 2);
  assert(one_plus_T.get_coeff(slong{0}).to<slong>() ==  1);
  assert(one_plus_T.get_coeff(slong{1}).is_zero());
  assert(one_plus_T.get_coeff(slong{2}).to<slong>() == -1);
}

int main() {
  test_default_ctor_is_zero();
  test_construction_and_raw();
  test_set_zero_one();
  test_arithmetic();
  test_modulo();
  test_comparison_and_map();
  test_divisibility();
  test_pow();
  test_to_conversion();
  test_ostream();
  test_zzx_basic();
  test_zzx_arithmetic();
  return 0;
}
