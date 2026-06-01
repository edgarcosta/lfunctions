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
  ZZ()                       { fmpz_init(z); }
  ZZ(slong n)                { fmpz_init(z); fmpz_set_si(z, n); }
  ZZ(ulong n)                { fmpz_init(z); fmpz_set_ui(z, n); }
  ZZ(const ZZ& o)            { fmpz_init(z); fmpz_set(z, o.z); }
  // Move: init the target to 0 first, then swap. This leaves the moved-from
  // object in a valid (zero) state that fmpz_clear can safely handle.
  ZZ(ZZ&& o) noexcept        { fmpz_init(z); fmpz_swap(z, o.z); }
  ~ZZ()                      { fmpz_clear(z); }

  ZZ& operator=(const ZZ& o)     { fmpz_set(z, o.z); return *this; }
  ZZ& operator=(ZZ&& o) noexcept { fmpz_swap(z, o.z); return *this; }
  ZZ& operator=(slong n)         { fmpz_set_si(z, n); return *this; }
  ZZ& operator=(ulong n)         { fmpz_set_ui(z, n); return *this; }

  fmpz*       _fmpz()       { return z; }
  const fmpz* _fmpz() const { return z; }

  bool is_zero() const { return fmpz_is_zero(z); }
};

} // namespace lfun
