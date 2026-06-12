# Public C API guide

`liblfun` exposes a single opaque-pointer C interface, declared in
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h).
This page is the narrative guide to that interface: the object lifecycle, the
mu / normalisation convention, how to supply Euler factors, how to read the
results back, the error model, and the hard limits. The complete
per-function reference, generated from the header, is at the
[bottom of this page](#full-header-reference); the prose below links into it
with cross-references such as {c:func}`Lfunc_init`.

Everything the library computes is returned as a *certified ball*: FLINT/Arb
`arb_t` (real) and `acb_t` (complex) values whose radius is a rigorous error
bound, never a bare floating-point number. To consume a result you read its
midpoint and radius with the usual Arb accessors, or compare it against a
reference ball with `arb_overlaps` / `acb_overlaps` / `arb_contains`.

```{contents}
:local:
:depth: 2
```

## At a glance

The whole interface follows one lifecycle. Create an L-function object, tell it
the Euler factors of the L-series, run the computation, query whichever
quantities you need, then free it.

```c
#include <flint/acb_poly.h>
#include "glfunc.h"

Lerror_t ecode = ERR_SUCCESS;     // accumulate errors here
double mus[2] = {0.0, 1.0};       // Gamma_R shifts (algebraic normalisation)

// 1. create: degree 2, conductor 37, weight-2 normalisation 0.5
Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ecode);
if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

// 2. supply the Euler factors (callback over every prime <= Lfunc_nmax(L))
ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

// 3. run the full computation
ecode |= Lfunc_compute(L);
if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

// 4. query
int64_t   rank  = Lfunc_rank(L);          // analytic rank
acb_srcptr sign = Lfunc_sign(L);          // root number epsilon
arb_srcptr lead = Lfunc_Taylor(L);        // leading Taylor coefficient
arb_srcptr zer  = Lfunc_zeros(L, 0);      // zeros of L on the critical line

// 5. free
Lfunc_clear(L);

// any warnings collected along the way:
fprint_errors(stderr, ecode);
```

The numbered steps map onto the sections below. The two helpers that bracket
the lifecycle, {c:func}`fatal_error` (does this accumulator hold a fatal flag?)
and {c:func}`fprint_errors` (print one human-readable line per flag), are
covered under [Error handling](#error-handling).

## The Lfunc_t lifecycle

### 1. Create: `Lfunc_init` and `Lfunc_init_advanced`

{c:func}`Lfunc_init` is the common entry point:

```c
Lfunc_t Lfunc_init(uint64_t degree, uint64_t conductor,
                   double normalisation, const double *mus,
                   Lerror_t *ecode);
```

- `degree` is the degree of the L-function. It must satisfy
  `2 <= degree <= MAX_DEGREE`; otherwise {c:macro}`ERR_BAD_DEGREE` fires.
- `conductor` is the (integer) conductor.
- `normalisation` and `mus` set the gamma factor and the algebraic-to-analytic
  shift; see [The mu / normalisation convention](#mu-normalisation).
  `mus` is an array of `degree` doubles that the library copies, so the caller
  may free it immediately after the call.
- `ecode` points at a {c:type}`Lerror_t` accumulator. On any fatal flag the
  returned pointer is unusable; check {c:func}`fatal_error` before continuing.

{c:func}`Lfunc_init_advanced` takes the same information through a
{c:struct}`Lparams_t` struct and additionally exposes the precision and
strategy knobs:

```c
typedef struct {
  uint64_t degree;
  uint64_t conductor;
  double   normalisation;
  double  *mus;
  int64_t  target_prec;   // target output precision, in bits
  int64_t  wprec;         // working precision, in bits
  int64_t  gprec;         // precision for the gamma-factor (G) computation
  int      self_dual;     // DK / YES / NO
  int      rank;          // DK, or a known rank
  char    *cache_dir;     // directory for the cached G data
  int      extract_powers; // YES to extract and assemble clean powers from exact fmpz supply
} Lparams_t;
```

`Lparams_t` is intentionally not a stable ABI or source-compatibility boundary.
Compile callers against the matching header and zero-initialize the struct before
setting fields; new fields may be added without preserving older struct layouts.

| Field | Meaning | Default / zero behavior |
| --- | --- | --- |
| {c:member}`target_prec <Lparams_t.target_prec>` | precision the results are refined to | {c:macro}`DEFAULT_TARGET_PREC` (100 bits) |
| {c:member}`wprec <Lparams_t.wprec>` | internal working precision | derived from `target_prec` |
| {c:member}`gprec <Lparams_t.gprec>` | precision of the gamma-factor product | derived from the working precision |
| {c:member}`self_dual <Lparams_t.self_dual>` | whether the L-function equals its dual | `DK`: the library decides |
| {c:member}`rank <Lparams_t.rank>` | analytic rank, if already known | `DK`: the library determines it |
| {c:member}`cache_dir <Lparams_t.cache_dir>` | where to read/write cached G data | current directory; see [The G-data cache](#g-data-cache) |
| {c:member}`extract_powers <Lparams_t.extract_powers>` | opt in to computing a certified pure power `L = M^k` by extracting `M` from exact `fmpz_poly_t` Euler-factor supply and assembling the result | `NO`: reject powers with {c:macro}`ERR_POWER` |

The tri-state fields use the macros {c:macro}`DK` (-1, "don't know"),
{c:macro}`YES` (1), and {c:macro}`NO` (0). Setting `self_dual = YES` when you
know the L-function is self-dual (every elliptic curve over Q, for instance)
lets the library skip computing the dual side. Supplying a known `rank` lets it
skip the rank search; if the computed rank then disagrees with what you
supplied, you get the {c:macro}`ERR_CONFLICT_RANK` warning.

The power/repeated-factor guard is enabled by default. If the supplied
L-function appears to be a perfect power, or to have repeated primitive factors,
the default behavior is to stop with the fatal {c:macro}`ERR_POWER` flag rather
than continue into a zero search with repeated zeros. Set
{c:member}`extract_powers <Lparams_t.extract_powers>` to {c:macro}`YES` only
when you want the library to certify a clean pure power `L = M^k`, compute the
primitive factor `M`, and assemble the rank, zeros, sign, Taylor coefficient,
and special values of `L` from it. The extracted factor is available through
{c:func}`Lfunc_factors`. This is not a general factorisation API: a primitive
object reports no factors, and an extracted pure power reports one primitive
factor with its multiplicity.

Extraction is certified only from exact integer-polynomial Euler factors supplied
through {c:func}`Lfunc_use_lpoly_fmpz` or
{c:func}`Lfunc_use_all_lpolys_fmpz`. Factors supplied through the `acb_poly_t`
APIs are still used for computation and for the repeated-factor guard, but they
do not carry exact provenance; if the guard detects a power from acb-only or
mixed acb/fmpz supply under `extract_powers = YES`, computation stops with
{c:macro}`ERR_POWER` rather than extracting.

Leaving a numeric knob at 0 asks the library to choose it. There is no need to
set `target_prec`, `wprec`, or `gprec` unless you have a specific reason;
{c:func}`Lfunc_init` is exactly `Lfunc_init_advanced` with those left to their
defaults, `self_dual = DK`, `rank = DK`, and no cache directory.

### 2. Supply the Euler factors

Between creation and computation you hand the library the local L-factors of
the L-series. The routes are described in detail in
[Supplying Euler factors](#supplying-euler-factors). In brief: either let the
library drive a callback over every prime it needs
({c:func}`Lfunc_use_all_lpolys` or
{c:func}`Lfunc_use_all_lpolys_fmpz`), or push factors one prime at a time
({c:func}`Lfunc_use_lpoly` or {c:func}`Lfunc_use_lpoly_fmpz`).

### 3. Compute

{c:func}`Lfunc_compute` runs the whole numerical pipeline (the FFT, the
upsampling, the zero search, the RH check, and the rank determination) and
returns an accumulated {c:type}`Lerror_t`. OR it into your accumulator and
check {c:func}`fatal_error` before reading any result. After it returns,
{c:func}`Lfunc_wprec` reports the working precision the computation actually
used.

### 4. Query

Read back the rank, root number, leading Taylor coefficient, zeros, special
values, and plot samples. See [Querying the results](#querying-the-results).

### 5. Clear

{c:func}`Lfunc_clear` frees the entire object. An {c:struct}`Lplot_t` returned
by {c:func}`Lfunc_plot_data` owns its own buffer and is freed separately with
{c:func}`Lfunc_clear_plot`.

(mu-normalisation)=

## The mu / normalisation convention

An L-function in this library is specified in the **algebraic** normalisation
(integer Dirichlet coefficients, functional equation `s -> k - s` for motivic
weight `k`) but is computed internally in the **analytic** normalisation
(functional equation centered at `1/2`). Two inputs bridge the two:

- `mus[i]` are the shifts of the `Gamma_R` factors that make up the completed
  L-function's gamma factor, in the algebraic normalisation.
- `normalisation` is the shift along the `s`-axis that converts algebraic to
  analytic. If `a_n` is the n-th Dirichlet coefficient in the algebraic
  normalisation, then `a_n / n^{normalisation}` is the coefficient in the
  analytic one. For a functional equation `Lambda(s) = eps Lambda(k - s)`, take
  `normalisation = (k - 1) / 2`.

Internally the library stores `mus[i] + normalisation` for each `i` and runs
entirely in the analytic normalisation. The hard constraint is:

> Every `mus[i] + normalisation` must be a **non-negative half-integer** (an
> element of `{0, 1/2, 1, 3/2, ...}`). Otherwise {c:func}`Lfunc_init` sets the
> fatal {c:macro}`ERR_MU_HALF` flag.

Because only the sum `mus[i] + normalisation` matters, the same L-function has
many equivalent `(mus, normalisation)` encodings. For an elliptic curve over Q
(degree 2, weight 2):

```
mus = [0, 1],   normalisation = 0.5     <==>     mus = [0.5, 1.5], normalisation = 0
```

Both feed the library the analytic shifts `[0.5, 1.5]` and must produce
identical output. For a classical modular form of weight 13 (motivic weight
12, so `normalisation = 6`):

```
mus = [0, 1],   normalisation = 6       <==>     mus = [6, 7],     normalisation = 0
```

Choosing `normalisation = (k - 1) / 2` and small `mus` keeps the encoding close
to how the object is usually tabulated (for example on the LMFDB), which is why
the worked examples use `mus = [0, 1]` with the weight-derived
`normalisation`.

## Supplying Euler factors

The library needs the local L-factor `L_p(T)` (a polynomial in `T` with the
constant term normalised to 1) for every prime `p` up to a bound fixed by the
conductor and degree. You may supply it numerically as an `acb_poly_t`, or
exactly as an integer `fmpz_poly_t`. {c:func}`Lfunc_nmax` returns that bound:

```c
uint64_t nmax = Lfunc_nmax(L);   // largest prime p for which a factor is expected
```

There are acb and exact-fmpz variants of the callback and one-prime-at-a-time
routes.

### Driven by the library: `Lfunc_use_all_lpolys`

```c
Lerror_t Lfunc_use_all_lpolys(
    Lfunc_t L,
    void (*lpoly_callback)(acb_poly_t lpoly, uint64_t p, int d,
                           int64_t prec, void *parm),
    void *param);
```

The library calls `lpoly_callback` once for each prime `p <= Lfunc_nmax(L)`,
in increasing order. The callback fills `lpoly` with `L_p(T)`; `d` is the
expected degree of the factor and `prec` the working precision, and `param` is
your opaque pointer passed straight through. This is the most convenient route:
you never enumerate the primes yourself.

It also carries a **short-circuit**: if the callback leaves `lpoly` set to the
zero polynomial for some prime, the library stops calling, takes that prime as
the end of your data, and reduces its internal `nmax` accordingly (as if you had
called {c:func}`Lfunc_reduce_nmax`). This is the idiomatic way to say "I have
factors up to here and no further": return zero for the first prime you cannot
supply.

A minimal callback that serves factors from a lookup table and zero-terminates
when it runs out:

```c
void lpoly_callback(acb_poly_t poly, uint64_t p, int d, int64_t prec, void *param)
{
  acb_poly_zero(poly);                 // default: zero -> short-circuits
  const factor_t *f = lookup(param, p);
  if (f)                               // we have a factor for this prime
    for (size_t i = 0; i < f->len; ++i)
      acb_poly_set_coeff_si(poly, i, f->coeff[i]);
}
```

### Exact integer callback: `Lfunc_use_all_lpolys_fmpz`

```c
Lerror_t Lfunc_use_all_lpolys_fmpz(
    Lfunc_t L,
    void (*lpoly_callback)(fmpz_poly_t lpoly, uint64_t p, int d, void *parm),
    void *param);
```

This is the same prime-driven interface, but the callback fills an exact
integer polynomial. Leaving `lpoly` equal to zero is the same insufficient-supply
sentinel. Use this route when you want `extract_powers = YES` to be able to
certify and assemble a pure power.

### One prime at a time: `Lfunc_use_lpoly`

```c
void Lfunc_use_lpoly(Lfunc_t L, uint64_t p, const acb_poly_t poly);
Lerror_t Lfunc_use_lpoly_fmpz(Lfunc_t L, uint64_t p, const fmpz_poly_t poly);
```

Push a single factor for the prime `p`. Use this when you already have the
primes enumerated (for example from `primesieve`) and want to loop yourself.
The library copies `poly`; the fmpz variant also converts it internally to acb
for the numerical computation while retaining the exact provenance needed for
power extraction. `Lfunc_use_lpoly_fmpz` returns any fatal retention error
immediately; callers should OR it into their accumulated `Lerror_t`. The
[worked example](#a-worked-end-to-end-example) on this page uses the callback form;
[`examples/rational.c`](https://github.com/edgarcosta/lfunctions/blob/main/examples/rational.c)
uses the exact push form over a `primesieve`-generated prime list.

### Supplying fewer factors than `nmax`: `Lfunc_reduce_nmax`

```c
bool Lfunc_reduce_nmax(Lfunc_t L, uint64_t nmax);
```

If you cannot reach {c:func}`Lfunc_nmax`, call this *before* supplying factors
to declare a smaller bound. **The library takes your word for it and does not
check** that the reduced bound still yields a correct result, so use it only
when you understand the consequence. Supplying too few factors (whether through
this call, or by zero-terminating the callback early, or by simply pushing
fewer with {c:func}`Lfunc_use_lpoly` / {c:func}`Lfunc_use_lpoly_fmpz`) is reported after the computation as the
{c:macro}`ERR_INSUFF_EULER` warning: the library expected more Euler factors
than it received, so the results may be degraded or, in the worst case,
meaningless.

```{note}
**Batch supply is on the way.** Array / batch supply of Euler factors and of
raw Dirichlet coefficients `a_n` (and the associated input-validation error
codes) is being added in the batch-supply API. Those entry points are not part
of the interface on `main` yet; this guide will document them once they land.
For now the supply routes are the callback and single-prime push APIs described
above.
```

## Querying the results

Call these only after {c:func}`Lfunc_compute` has returned without a fatal
flag. The `arb_srcptr` / `acb_srcptr` return values are **borrowed**: they point
into the `Lfunc_t` and stay valid until {c:func}`Lfunc_clear`. Do not free
them; copy out (`arb_set` / `acb_set`) anything you need to outlive the object.

### Rank: `Lfunc_rank`

{c:func}`Lfunc_rank` returns the analytic rank (order of vanishing at the
central point) as an `int64_t`. A returned rank of 0 or 1 is rigorous; a value
greater than 1 is not certified. If the rank could not be determined you get
the {c:macro}`ERR_NO_RANK` warning, and if it disagrees with a rank you supplied
through {c:struct}`Lparams_t` you get {c:macro}`ERR_CONFLICT_RANK`.

### Root number: `Lfunc_sign` and `Lfunc_sqrt_sign`

{c:func}`Lfunc_sign` returns the sign (root number) `eps` of the functional
equation `Lambda(s) = eps * Lambda(k - s)`, as a borrowed `acb_srcptr`. It is a
complex number of absolute value 1; `|eps| = 1` is a useful sanity check on a
computation.

```{note}
The relevant symbol is {c:func}`Lfunc_sign`. Some example programs print this
quantity under the label `Epsilon`, but there is no `Lfunc_epsilon` function;
the accessor is `Lfunc_sign`.
```

{c:func}`Lfunc_sqrt_sign` returns a square root of the sign, with the branch
chosen so that `Lambda` is positive at the central point.

### Leading Taylor coefficient: `Lfunc_Taylor`

{c:func}`Lfunc_Taylor` returns the first non-zero Taylor coefficient at the
central point, that is `L^{(rank)}((w + 1)/2) / rank!`, where `w/2` is the
normalisation (when the input is algebraic, `w` is the motivic weight). It
carries the same rigour caveat as the rank: rigorous for rank 0 or 1. To
evaluate the L-function exactly at the central point, use this rather than
{c:func}`Lfunc_special_value`.

### Extracted factors: `Lfunc_factors`

```c
uint64_t Lfunc_factors(Lfunc_t L, Lfunc_t **factors, uint64_t **mults);
```

For a primitive object, or for an object that was not produced through
{c:member}`extract_powers <Lparams_t.extract_powers>`, this returns 0. For a
successfully extracted pure power `L = M^k`, it returns 1. If `factors` is
non-NULL, `*factors` is set to a borrowed array whose first element points to
`M`; if `mults` is non-NULL, `*mults` is set to an array whose first element is
`k`. The returned arrays and factor objects are owned by `L`; do not clear them
yourself, and do not use them after {c:func}`Lfunc_clear(L)`.

The assembled object exposes the quantities that can be assembled from the
factor: rank, sign, zeros, leading Taylor coefficient, and special values.
Plot samples are the exception; see [Plot data](#plot-data-lfunc-plot-data-and-lfunc-clear-plot).

### Zeros: `Lfunc_zeros`

```c
arb_srcptr Lfunc_zeros(Lfunc_t L, uint64_t side);
```

Returns the array of imaginary parts of the zeros on the critical line, as
borrowed `arb_t` balls. `side = 0` gives the zeros of `L`; `side = 1` gives
those of the dual L-function (for a self-dual L-function the two agree). When
the rank is 0 or 1 and the RH check succeeded (no {c:macro}`ERR_RH_ERROR`
warning), the list is complete up to the height reached; otherwise some zeros
may be missing. At most {c:macro}`MAX_ZEROS` (256) zeros are stored per side.
If the zeros could not be refined to the target precision you get the
{c:macro}`ERR_ZERO_PREC` warning.

### Special values: `Lfunc_special_value`

```c
Lerror_t Lfunc_special_value(acb_t res, Lfunc_t L, double re, double im);
```

Writes `L(re + i*im)` into `res` (an `acb_t` the caller has initialised) and
returns an accumulated {c:type}`Lerror_t`. Accuracy is excellent near the
critical line and falls away rapidly as you move away from it, but values such
as `L(k)` and `L(0)` come out sensibly. The routine requires `im >= 0`;
`im < 0` raises {c:macro}`ERR_SPEC_NZ`. For the central point use
{c:func}`Lfunc_Taylor` instead. The lower-level
{c:func}`Lfunc_special_value_choice` additionally returns the derivative and
lets you choose between `L` and the completed `Lambda`. For an assembled power,
`do_dash = true` is supported only for the completed function (`lam_p = true`);
the non-completed derivative returns {c:macro}`ERR_SPEC_VALUE`.

### Plot data: `Lfunc_plot_data` and `Lfunc_clear_plot`

```c
Lplot_t *Lfunc_plot_data(Lfunc_t L, uint64_t side, double max_t,
                         uint64_t n_points);
void     Lfunc_clear_plot(Lplot_t *Lp);
```

Returns roughly `n_points` samples of the real-valued Z-function
`exp(i*theta(t)) * L(k/2 + i*t)` over `t` in `[0, max_t]`, for `L` (`side = 0`)
or its dual (`side = 1`). The samples come back as plain doubles in a
heap-allocated {c:struct}`Lplot_t`:

```c
typedef struct {
  uint64_t n_points;   // number of samples actually produced
  double   spacing;    // t-step between consecutive samples
  double  *points;     // the sample values, length n_points
} Lplot_t;
```

Free it with {c:func}`Lfunc_clear_plot` when done. Unlike the query accessors
above, this struct is a separate allocation that you own.

```{warning}
Plot samples are not meaningful for an object assembled through
{c:member}`extract_powers <Lparams_t.extract_powers>`, because that path skips
the sample-generating pipeline. Use {c:func}`Lfunc_factors` to obtain the
primitive factor and plot that factor instead.
```

## Error handling

Every API entry point that can fail returns an {c:type}`Lerror_t`, a 64-bit
bitfield. The convention is:

- **The lower 32 bits are fatal**; if any is set the computation cannot
  continue and the results are invalid.
- **The upper 32 bits are warnings**; the computation produced output, but with
  a caveat (degraded precision, an unverified rank, an incomplete zero list).

The idiom is to keep one accumulator and OR every return value into it, then
gate on {c:func}`fatal_error` before each next step:

```c
Lerror_t ecode = ERR_SUCCESS;
ecode |= Lfunc_use_all_lpolys(L, cb, NULL);
ecode |= Lfunc_compute(L);
if (fatal_error(ecode)) {       // true iff any low-32 (fatal) bit is set
  fprint_errors(stderr, ecode); // one human-readable line per set flag
  /* abort this L-function */
}
/* ... query ... */
fprint_errors(stderr, ecode);   // surface any warnings collected en route
```

{c:func}`fatal_error` returns `true` iff any fatal (low-32) bit is set;
warnings alone return `false`. {c:func}`fprint_errors` writes one line per set
flag to the given `FILE *`. {c:macro}`ERR_SUCCESS` (0) is the clean state.

### Fatal codes

| Macro | Meaning |
| --- | --- |
| {c:macro}`ERR_NO_DATA` | the first two `Lambda(t)` samples contained zero: nothing usable |
| {c:macro}`ERR_ZERO_ERROR` | unexpected failure isolating zeros |
| {c:macro}`ERR_OOM` | out of memory |
| {c:macro}`ERR_UPSAMPLE` | failure during upsampling |
| {c:macro}`ERR_MU_HALF` | some `mus[i] + normalisation` is not a non-negative half-integer |
| {c:macro}`ERR_M_ERROR` | could not compute the coefficient-count bound `M` |
| {c:macro}`ERR_STAT_POINT` | failure locating a stationary point |
| {c:macro}`ERR_SPEC_VALUE` | failure in the special-value routine |
| {c:macro}`ERR_G_INFILE` | failure reading the G-data cache file |
| {c:macro}`ERR_BAD_DEGREE` | degree below 2 or above {c:macro}`MAX_DEGREE` |
| {c:macro}`ERR_SPEC_NZ` | special value requested with `Im(s) < 0` |
| {c:macro}`ERR_G_EXTENT` | the G grid does not extend low enough (conductor too large for the fixed grid floor, or a stale cached grid was reused) |
| {c:macro}`ERR_POWER` | the supplied L-function is a perfect power or has repeated primitive factors; set {c:member}`extract_powers <Lparams_t.extract_powers>` to {c:macro}`YES` and supply exact fmpz Euler factors to request certified extraction of a clean pure power |

### Warning codes

| Macro | Meaning |
| --- | --- |
| {c:macro}`ERR_SOME_DATA` | had sensible data, but not all the way through the Turing zone |
| {c:macro}`ERR_ZERO_PREC` | could not isolate zeros to the target precision |
| {c:macro}`ERR_RH_ERROR` | the RH check on the computed zeros failed |
| {c:macro}`ERR_INSUFF_EULER` | ran out of Euler factors before the expected bound |
| {c:macro}`ERR_NO_RANK` | could not determine the rank |
| {c:macro}`ERR_CONFLICT_RANK` | the computed rank disagreed with the supplied one |
| {c:macro}`ERR_DBL_ZERO` | a stationary point failed to converge (a double zero?) |
| {c:macro}`ERR_SPEC_PREC` | could not reach the target error bound in a special value |
| {c:macro}`ERR_G_OUTFILE` | could not open the file to write the G-data cache |

When you add a new error code, it belongs in the header next to these and needs
a matching message branch in `fprint_errors` (in `src/glfunc.c`).

## Hard limits

```{list-table}
:header-rows: 1

* - Macro
  - Value
  - Meaning
* - {c:macro}`MAX_DEGREE`
  - 9
  - largest supported degree; raising it requires extending the Buthe integral
    tables in `src/buthe.c` and reviewing `src/g.c`
* - {c:macro}`MAX_ZEROS`
  - 256
  - hard cap on the number of zeros found and stored per side
* - {c:macro}`DEFAULT_TARGET_PREC`
  - 100
  - default target precision, in bits, when {c:member}`target_prec <Lparams_t.target_prec>` is left at 0
```

({c:macro}`MAX_R` is an alias for {c:macro}`MAX_DEGREE` that sizes the Buthe
tables.)

(g-data-cache)=

## The G-data cache and `cache_dir`

The gamma-factor product (the "G data") depends only on the degree, the gamma
shifts, and the precision, not on the Euler factors. Computing it is expensive,
so the library caches it to disk, keyed by the analytic normalisation, in a file
named `g_<normalisation>` (for example `g_0.5`). On a subsequent run with a
matching gamma factor the cache is read back instead of recomputed.

{c:member}`cache_dir <Lparams_t.cache_dir>` (settable through
{c:func}`Lfunc_init_advanced`) chooses the directory for these files; left
unset, the current working directory is used.

```{warning}
A stale or mismatched `g_<normalisation>` file in the cache directory can poison
a run. For a hermetic computation, point {c:member}`cache_dir
<Lparams_t.cache_dir>` at a clean directory, or remove any `g_*` files from the
working directory first. A cache that cannot be read raises the fatal
{c:macro}`ERR_G_INFILE`; a grid that does not extend far enough (which a stale
cache can also cause) raises {c:macro}`ERR_G_EXTENT`.
```

## A worked end-to-end example

The example below builds the degree-2 L-function of the Ramanujan tau form
(the unique newform of weight 12 and level 1) from its Euler factors and reads
back the rank, root number, leading Taylor coefficient, two special values, and
the first zeros. It is a condensed version of
[`examples/tau.cpp`](https://github.com/edgarcosta/lfunctions/blob/main/examples/tau.cpp),
whose committed expected output the values below are taken from, so they are
verified rather than illustrative. The motivic weight is 11, hence
`normalisation = 5.5` and `mus = [0, 1]`.

```c
#include <flint/acb_poly.h>
#include <stdio.h>
#include "glfunc.h"

// A few good Euler factors L_p(T) = 1 - a_p T + p^{11} T^2, from the LMFDB
// page for the weight-12 level-1 newform 1.12.a.a. (The real program carries
// every prime up to Lfunc_nmax; this excerpt keeps the table short.)
static const struct { uint64_t p; long c[3]; } table[] = {
  { 2, {1,      24,             2048} },
  { 3, {1,    -252,           177147} },
  { 5, {1,   -4830,         48828125} },
  { 7, {1,   16744,       1977326743} },
  {11, {1, -534612,     285311670611} },
  /* ... up to Lfunc_nmax(L) ... */
};

static void lpoly_callback(acb_poly_t poly, uint64_t p, int d,
                           int64_t prec, void *param)
{
  acb_poly_zero(poly);                 // unknown prime -> zero -> short-circuit
  for (size_t i = 0; i < sizeof(table)/sizeof(table[0]); ++i)
    if (table[i].p == p) {
      for (int j = 0; j < 3; ++j)
        acb_poly_set_coeff_si(poly, j, table[i].c[j]);
      return;
    }
}

int main(void)
{
  Lerror_t ecode = ERR_SUCCESS;
  double   mus[2] = {0.0, 1.0};

  // degree 2, conductor 1, motivic weight 11 -> normalisation 5.5
  Lfunc_t L = Lfunc_init(2, 1, 5.5, mus, &ecode);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

  ecode |= Lfunc_compute(L);
  if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

  printf("rank = %ld\n", (long) Lfunc_rank(L));   // 0

  printf("epsilon = ");                            // 1 (this form is self-dual)
  acb_printd(Lfunc_sign(L), 20); printf("\n");

  printf("leading Taylor coeff = ");               // 0.79212283864603056936...
  arb_printd(Lfunc_Taylor(L), 20); printf("\n");

  acb_t v; acb_init(v);
  ecode |= Lfunc_special_value(v, L, 6.5, 0.0);    // 0.83934551203194208649...
  printf("L(6.5) = "); acb_printd(v, 20); printf("\n");
  acb_clear(v);

  arb_srcptr zeros = Lfunc_zeros(L, 0);            // first zero 9.2223793999...
  printf("first zero = "); arb_printd(zeros + 0, 20); printf("\n");

  Lfunc_clear(L);
  fprint_errors(stderr, ecode);                    // print any warnings
  return 0;
}
```

Running the full example prints, among other lines:

```text
Order of vanishing = 0
Epsilon = (1 + 0j)  +/-  (1.09e-105, 4.67e-53j)
First non-zero Taylor coeff = 0.79212283864603056936 +/- 6.0749e-52
L(6.5) = (0.83934551203194208649 - 2.1002674689610803795e-51j)  +/-  (1.29e-37, 1.29e-37j)
Zero 0 = 9.2223793999211025222 +/- 2.1895e-47
```

Each printed value is a ball: the `+/-` part is the rigorous radius. A test
asserts correctness with `acb_overlaps` / `arb_overlaps` against reference
balls rather than by string-matching this output.

### More examples in the tree

The [`examples/`](https://github.com/edgarcosta/lfunctions/tree/main/examples)
directory has runnable programs covering the API end to end. A few worth
reading:

- [`tau.cpp`](https://github.com/edgarcosta/lfunctions/blob/main/examples/tau.cpp):
  the full version of the example above (Ramanujan tau, degree 2).
- [`ec_37.a1.cpp`](https://github.com/edgarcosta/lfunctions/blob/main/examples/ec_37.a1.cpp):
  a rank-1 elliptic curve, also exercising {c:func}`Lfunc_sqrt_sign` and a
  BSD leading-coefficient check.
- [`rational.c`](https://github.com/edgarcosta/lfunctions/blob/main/examples/rational.c):
  a generic command-line driver that parses
  `label:degree:conductor:weight:[mus]:[[euler_factors]]` lines and uses the
  {c:func}`Lfunc_use_lpoly_fmpz` push route over a `primesieve` prime list.

## Full header reference

The declarations below are generated from
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h)
by Doxygen and rendered through Breathe, so the header stays the single source
of truth for the per-function reference.

```{doxygenfile} glfunc.h
:project: liblfun
```
