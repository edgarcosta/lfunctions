// Unit tests for lfun::ZZ and lfun::ZZX (lflint.h).
// Plain <cassert>; exits 0 on success.
#include <cassert>
#include "lflint.h"

using lfun::ZZ;

static void test_default_ctor_is_zero() {
  ZZ a;
  assert(a.is_zero());
}

int main() {
  test_default_ctor_is_zero();
  return 0;
}
