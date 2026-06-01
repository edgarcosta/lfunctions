// Unit tests for lfun::ZZ and lfun::ZZX (lflint.h).
// Plain <cassert>; exits 0 on success.
#include <cassert>
#include <utility>
#include "lflint.h"

using lfun::ZZ;

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
}

int main() {
  test_default_ctor_is_zero();
  test_construction_and_raw();
  test_set_zero_one();
  test_arithmetic();
  return 0;
}
