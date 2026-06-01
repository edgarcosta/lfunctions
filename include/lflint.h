// Copyright Edgar Costa 2026
// See LICENSE file for license details.
//
// lflint.h — minimal RAII C++ wrappers over the FLINT C API.
//
// Two value types in namespace lfun:
//   ZZ   — wraps fmpz_t       (arbitrary-precision integer)
//   ZZX  — wraps fmpz_poly_t  (integer polynomial in T)
//
// Covers exactly the surface used by this repo's examples; not a
// general-purpose flintxx replacement. See
// docs/superpowers/specs/2026-05-29-drop-flintxx-design.md.

#pragma once

#include <ostream>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

namespace lfun {

class ZZ {
  fmpz_t z;
public:
  ZZ() { fmpz_init(z); }
  ~ZZ() { fmpz_clear(z); }

  bool is_zero() const { return fmpz_is_zero(z); }
};

} // namespace lfun
