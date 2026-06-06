# Symmetric-power L-functions of an elliptic curve

This is a how-to for computing the symmetric-power L-function `Sym^k(E)` of an
elliptic curve `E/Q` with `lfunctions`. The good-prime Euler factors are formed
from the curve's `a_p` by the helper `sym_power_lpoly`
([`include/sym_power.h`](../include/sym_power.h)); a complete worked example
(Sym^2 and Sym^3 of `11.a1`, with golden-value regression assertions) is
[`examples/ec_sym.cpp`](../examples/ec_sym.cpp).

There is **no native `Sym^k` constructor yet** (tracked as a future feature):
you assemble the L-function by hand from its degree, normalisation, `mus`,
`self_dual` flag, and Euler factors, exactly as for any object built through the
public API. For the general C API (init, supplying factors, computing, querying)
see [`doc/api.md`](api.md); this document covers only the `Sym^k`-specific
conventions and the factor assembly.

## Conventions

For `Sym^k` of an elliptic curve over Q (motivic weight `k`):

| Parameter | Value | Notes |
| --- | --- | --- |
| `degree` | `k + 1` | `Sym^1` is the curve itself (degree 2). |
| `normalisation` | `k / 2` | The algebraic-to-analytic shift; `Lfunc_init`'s `normalisation` argument, or `Lparams_t.normalisation`. |
| `self_dual` | `YES` | A symmetric power of the self-dual L-function of an elliptic curve is self-dual. |
| `mus` | length `k + 1`, see below | The Gamma_R shifts. |

The `mus` vector has `k + 1` entries. The closed form (the one the high-degree
suite uses) is:

```
u = ceil(k / 2)
mus[i]     = -i        for i in [0, u)
mus[i + u] = -i + 1    for i in [0, u)
if k is even (so degree = k + 1 > k):  mus[k] = -2 * floor(u / 2)
```

Worked out for the first few powers:

| `k` | `degree` | `normalisation` | `mus` |
| --- | --- | --- | --- |
| 1 | 2 | 0.5 | `[0, 1]` (the curve itself) |
| 2 | 3 | 1.0 | `[0, 1, 0]` |
| 3 | 4 | 1.5 | `[0, -1, 1, 0]` |
| 4 | 5 | 2.0 | `[0, -1, 1, 0, -2]` |

(`mus[i] + normalisation` must be a non-negative half-integer, which these
satisfy.) The `k = 2` and `k = 3` rows are exactly the vectors used in
`examples/ec_sym.cpp`.

## Forming the good-prime Euler factors

At a prime `p` of **good** reduction the Satake parameters `alpha, beta` of `E`
satisfy `alpha + beta = a_p` and `alpha * beta = p`. The local factor of the
k-th symmetric power is the integer polynomial of degree `k + 1`

```
out(T) = prod_{i=0}^{k} (1 - alpha^i beta^{k-i} T).
```

`sym_power_lpoly` computes this exactly, without complex roots, from the Lucas
sequence `V_m = alpha^m + beta^m` (`V_0 = 2`, `V_1 = a_p`,
`V_m = a_p*V_{m-1} - p*V_{m-2}`):

```c
#include "sym_power.h"   // include/sym_power.h

void sym_power_lpoly(fmpz_poly_t out, slong a_p, ulong p, int k);
```

- **Input:** the curve's `a_p`, the good prime `p`, and the power `k >= 1`.
- **Output:** `out`, the degree-`(k + 1)` local L-polynomial in ascending order
  (constant term first), with constant term 1.
- With `k = 1` this returns the curve's own factor `1 - a_p*T + p*T^2`.

The output is an `fmpz_poly_t` (arbitrary-precision integer coefficients) on
purpose: the coefficients routinely exceed 64 bits (the leading one has absolute
value `p^{k(k+1)/2}`), so an `int64` representation would overflow from `Sym^5`
upward. Convert it to the `acb_poly_t` the library expects with
`acb_poly_set_fmpz_poly(poly, out, prec)`.

### a_p are caller-supplied

The library has no elliptic-curve point counting and no Hecke machinery, so
**you must supply the `a_p` yourself** (from the LMFDB, smalljac/`lpdata`,
Pari/GP, Sage, or Magma). `sym_power_lpoly` is a pure function of `(a_p, p, k)`;
it turns your `a_p` into the `Sym^k` factor but does not compute `a_p`. You need
`a_p` at every good prime up to `Lfunc_nmax` (query it after init), which grows
with `k` because `nmax` scales with the conductor of `Sym^k(E)`.

A convenient pattern is to drive the supply step with `Lfunc_use_all_lpolys`,
whose callback is invoked for each prime `p <= Lfunc_nmax`: look up `a_p`, call
`sym_power_lpoly`, convert to `acb_poly_t`, and return it. (At a missing good
prime, returning the zero polynomial only OR-s in the warning
`ERR_INSUFF_EULER`, not a fatal error, so a gap silently corrupts the result;
`examples/ec_sym.cpp` deliberately aborts loudly instead.)

## The bad-factor caveat

`sym_power_lpoly` is for **good primes only**. Bad primes are the caller's
responsibility, and they split into two cases:

- **Multiplicative reduction.** At a prime of multiplicative reduction the
  `Sym^k` bad factor is the degree-1 polynomial `1 - a_p^k * T`, where
  `a_p = +1` (split) or `a_p = -1` (non-split). This **is** a function of `a_p`,
  so you can write it down directly. For example `11.a1` has split
  multiplicative reduction at 11 with `a_11 = +1`, so its `Sym^k` factor at 11 is
  `1 - T` (coefficients `[1, -1]`) for every `k`. Pad the factor to `degree + 1`
  coefficients with trailing zeros if your supply path expects a fixed width.
- **Additive (ramified) reduction.** At a prime of additive reduction the
  `Sym^k` bad factor is **not** a function of `a_p` and cannot be derived from
  it. (For instance, `Sym^4` of `27.a` at `p = 3` is `1 - p^2*T`, which is not
  recoverable from `a_3 = 0`.) These factors are out of scope for the helper and
  must be supplied from a reliable source such as the LMFDB.

`examples/ec_sym.cpp` deliberately uses `11.a1`, whose only bad prime is
multiplicative, to stay within the part of the workflow the helper supports.

## Putting it together

The end-to-end recipe, as carried out in `examples/ec_sym.cpp`:

1. Choose `k`. Set `degree = k + 1`, `normalisation = k / 2`, `self_dual = YES`,
   and the `mus` vector above; initialise with `Lfunc_init_advanced`.
2. Query `Lfunc_nmax`. Make sure you have `a_p` at every good prime up to it.
3. Supply the factors (for example via `Lfunc_use_all_lpolys`): at each good
   prime, `sym_power_lpoly(out, a_p, p, k)` then
   `acb_poly_set_fmpz_poly`; at each bad prime, the explicit bad factor above.
4. `Lfunc_compute`, then query rank, root number, zeros, leading Taylor
   coefficient, and special values as usual (see [`doc/api.md`](api.md)).

For degree `>= 3` (so any `Sym^k` with `k >= 2`) the computation currently emits
a tolerated `ERR_RH_ERROR` **warning** (the Turing zero count is miscalibrated
for degree `>= 3`); it is a warning, not a fatal error, so collect and report it
rather than asserting against it, as `examples/ec_sym.cpp` does. While that
warning is present, rank/zero checks in this workflow are comparisons against
trusted golden data, not RH-certified rank/zero-completeness claims.

The same `sym_power_lpoly` implementation backs both `examples/ec_sym.cpp` and
the high-degree regression driver, and is checked against Pari `lfunsympow` for
`k = 1..8` in `test/sympow_lpoly_test.c`.
