# Public C API guide

`liblfun` exposes a single opaque-pointer C interface, declared in
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h).
This page is the narrative guide to that interface: the object lifecycle, the
mu / normalisation convention, how to supply Euler factors or Dirichlet
coefficients, how to read the results back, the error model, and the hard
limits. The complete
per-function reference, generated from the header, is at the
[bottom of this page](#full-header-reference); the prose below links into it
with cross-references such as {c:func}`Lfunc_init`.

Primary numerical query results are returned as *certified balls*: FLINT/Arb
`arb_t` (real) and `acb_t` (complex) values whose radius is a rigorous error
bound. Plot samples are the exception: {c:func}`Lfunc_plot_data` returns plain
doubles for visualisation. To consume a ball result you read its midpoint and
radius with the usual Arb accessors, or compare it against a reference ball with
`arb_overlaps` / `acb_overlaps` / `arb_contains`.

```{contents}
:local:
:depth: 2
```

## At a glance

The whole interface follows one lifecycle. Create an L-function object, supply
either local Euler factors or raw Dirichlet coefficients for the L-series, run
the computation, query whichever quantities you need, then free it.

```c
#include <flint/acb_poly.h>
#include "glfunc.h"

Lerror_t ecode = ERR_SUCCESS;     // accumulate errors here
double mus[2] = {0.0, 1.0};       // Gamma_R shifts (algebraic normalisation)

// 1. create: degree 2, conductor 37, weight-2 normalisation 0.5
Lfunc_t L = Lfunc_init(2, 37, 0.5, mus, &ecode);
if (fatal_error(ecode)) { fprint_errors(stderr, ecode); return 1; }

// 2. supply the Euler factors (callback over every prime <= Lfunc_nmax(L));
//    array supply and raw Dirichlet-coefficient supply are also available.
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
} Lparams_t;
```

{c:func}`Lfunc_init_advanced` copies most of these fields verbatim: it does
**not** re-apply {c:func}`Lfunc_init`'s defaults, so a zeroed
{c:struct}`Lparams_t` is not equivalent to {c:func}`Lfunc_init`. Only `wprec`
and `gprec` treat 0 as "derive a sensible value"; the other fields take their
value literally.

| Field | Meaning | Effect of the value you pass |
| --- | --- | --- |
| {c:member}`target_prec <Lparams_t.target_prec>` | precision the results are refined to | copied verbatim; pass {c:macro}`DEFAULT_TARGET_PREC` (100 bits) for the standard target. A literal 0 is a 0-bit target, not the default |
| {c:member}`wprec <Lparams_t.wprec>` | internal working precision | 0 means derive it from `target_prec` (the one field with a 0-default) |
| {c:member}`gprec <Lparams_t.gprec>` | precision of the gamma-factor product | 0 means derive it from the target/working precision |
| {c:member}`self_dual <Lparams_t.self_dual>` | whether the L-function equals its dual | copied verbatim: `DK` lets the library decide, `YES` skips the dual side, and the literal 0 is `NO` |
| {c:member}`rank <Lparams_t.rank>` | analytic rank, if already known | copied verbatim: `DK` lets the library determine it; a literal 0 asserts rank 0 (which can raise {c:macro}`ERR_CONFLICT_RANK`) |
| {c:member}`cache_dir <Lparams_t.cache_dir>` | where to read/write cached G data | copied verbatim; `NULL` disables the on-disk cache. See [The G-data cache](#g-data-cache) |

The tri-state fields use the macros {c:macro}`DK` (-1, "don't know"),
{c:macro}`YES` (1), and {c:macro}`NO` (0). Setting `self_dual = YES` when you
know the L-function is self-dual (every elliptic curve over Q, for instance)
lets the library skip computing the dual side. Supplying a known `rank` lets it
skip the rank search; if the computed rank then disagrees with what you
supplied, you get the {c:macro}`ERR_CONFLICT_RANK` warning.

Because the copy is verbatim, set `self_dual` and `rank` to {c:macro}`DK`
unless you really mean `NO` / rank 0, and pass
{c:macro}`DEFAULT_TARGET_PREC` for `target_prec` unless you have a specific
reason to change it; you can still leave `wprec` and `gprec` at 0 to have them
derived. {c:func}`Lfunc_init` is exactly `Lfunc_init_advanced` with
`target_prec = DEFAULT_TARGET_PREC`, `wprec = gprec = 0`,
`self_dual = DK`, `rank = DK`, and `cache_dir = "."`.

### 2. Supply factors or coefficients

Between creation and computation you hand the library the L-series data. The
routes are described in detail in
[Supplying factors and coefficients](#supplying-euler-factors). In brief: you
can let the library drive a callback over every prime it needs
({c:func}`Lfunc_use_all_lpolys`), push factors one prime at a time
({c:func}`Lfunc_use_lpoly`), pass a consecutive-prime array of factors
({c:func}`Lfunc_use_lpolys_fmpz` / {c:func}`Lfunc_use_lpolys_acb`), or pass raw
Dirichlet coefficients ({c:func}`Lfunc_use_dirichlet_coeffs_fmpz` /
{c:func}`Lfunc_use_dirichlet_coeffs_acb`).

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
(integer Dirichlet coefficients, functional equation `s -> w + 1 - s` for
motivic weight `w`) but is computed internally in the **analytic**
normalisation (functional equation centered at `1/2`). Two inputs bridge the
two:

- `mus[i]` are the shifts of the `Gamma_R` factors that make up the completed
  L-function's gamma factor, in the algebraic normalisation.
- `normalisation` is the shift along the `s`-axis that converts algebraic to
  analytic. If `a_n` is the n-th Dirichlet coefficient in the algebraic
  normalisation, then `a_n / n^{normalisation}` is the coefficient in the
  analytic one. For motivic weight `w`, take `normalisation = w / 2`. For a
  classical modular form of weight `k`, this is `(k - 1) / 2`.

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

Choosing the weight-derived `normalisation` and small `mus` keeps the encoding
close to how the object is usually tabulated (for example on the LMFDB), which
is why the worked examples use `mus = [0, 1]` with that convention.

(supplying-euler-factors)=

## Supplying factors and coefficients

The library needs enough L-series data up to a conductor- and degree-dependent
bound. For Euler-factor routes, the bound is the largest prime `p` for which a
local factor is expected; for raw coefficient routes, it is the largest
Dirichlet coefficient index `n` expected. {c:func}`Lfunc_nmax` returns that
bound:

```c
uint64_t nmax = Lfunc_nmax(L);
```

Euler-factor supply takes local L-polynomials `L_p(T)` (polynomials in `T`
with constant term normalised to 1) in the **algebraic** normalisation. Raw
Dirichlet-coefficient supply can take either algebraic or analytic
coefficients; the caller must say which one with {c:macro}`LFUNC_ALGEBRAIC_NORM` or
{c:macro}`LFUNC_ANALYTIC_NORM`.

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
your opaque pointer passed straight through. A nonzero factor must have constant
term 1 and degree at most the L-function degree. A `NULL` callback is fatal
{c:macro}`ERR_BAD_SUPPLY`. This is the most convenient route: you never
enumerate the primes yourself.

It also carries a **short-circuit**: if the callback leaves `lpoly` set to the
zero polynomial for some prime, the library stops calling, takes that prime as
the end of your data, reduces its internal `nmax` to `p - 1`, and returns the
warning bit {c:macro}`ERR_INSUFF_EULER`. This is the idiomatic way to say "I
have factors up to here and no further": return zero for the first prime you
cannot supply.

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

### One prime at a time: `Lfunc_use_lpoly`

```c
void Lfunc_use_lpoly(Lfunc_t L, uint64_t p, const acb_poly_t poly);
```

Push a single factor for the prime `p`. The prime must be at least 2, and
successive pushes must be in strictly increasing prime order. A non-prime `p`
is fatal {c:macro}`ERR_BAD_SUPPLY`. Use this when you
already have the primes enumerated (for example from `primesieve`) and want to
loop yourself. The factor must be nonzero, have exact constant term 1, and have
degree at most the L-function degree. If {c:func}`Lfunc_nmax` has not yet been called,
the first successful push computes and freezes that bound before multiplying the
factor into the coefficient array.

This route has no automatic end-of-data signal: if your explicit loop will stop
before the original bound, call {c:func}`Lfunc_reduce_nmax` before
{c:func}`Lfunc_compute`. Otherwise the library cannot distinguish "done" from
"more pushes are coming".

Because this function returns `void`, fatal supply errors such as `p < 2`,
non-prime `p`, a `NULL` or malformed `poly` ({c:macro}`ERR_BAD_SUPPLY`), a
duplicate/out-of-order push, or an attempt to push after another supply route
({c:macro}`ERR_SUPPLY_CONFLICT`) are recorded on the object and surfaced by
{c:func}`Lfunc_compute`.

The [worked example](#a-worked-end-to-end-example) on this page uses the
callback form; other examples use the push route when they already own an
explicit prime loop.

### Consecutive-prime arrays: `Lfunc_use_lpolys_fmpz` and `Lfunc_use_lpolys_acb`

```c
Lerror_t Lfunc_use_lpolys_fmpz(Lfunc_t L,
                               const fmpz_poly_struct *f,
                               uint64_t len);
Lerror_t Lfunc_use_lpolys_acb(Lfunc_t L,
                              const acb_poly_struct *f,
                              uint64_t len);
```

These batch front-ends supply one Euler factor per consecutive prime:
`f[0]` is the factor at `p = 2`, `f[1]` at `p = 3`, and so on. The caller does
not pass primes; the library sieves them and consumes the k-th supplied factor
for the k-th prime up to {c:func}`Lfunc_nmax`.

The factor arrays are in the same algebraic normalisation as
{c:func}`Lfunc_use_lpoly`, and internally route through the same multiplication
path. The `fmpz` form is for exact integer polynomials and is converted to
`acb_poly_t` at working precision; the `acb` form is the ball-valued route for
callers that already have certified polynomial balls.

Each supplied factor must be nonzero, have exact constant term 1, and have
degree at most the L-function degree. If `len` is shorter than the number of primes up to
`nmax`, the library stops at the first missing prime, reduces `nmax` to `p - 1`,
and returns the warning bit {c:macro}`ERR_INSUFF_EULER`; surplus factors are
ignored. A `NULL` factor array is allowed only when `len == 0`; with positive
`len` it is fatal {c:macro}`ERR_BAD_SUPPLY`.

Choose exactly one Euler-factor route for an object: callback, push, or one
factor-array call. Repeated factor-array calls and route mixing are fatal
{c:macro}`ERR_SUPPLY_CONFLICT`. They are not append/chunk operations, and
separate same-prime factors cannot be multiplied after the local series have
already been expanded.

### Raw coefficient arrays: `Lfunc_use_dirichlet_coeffs_fmpz` and `Lfunc_use_dirichlet_coeffs_acb`

```c
Lerror_t Lfunc_use_dirichlet_coeffs_fmpz(Lfunc_t L,
                                         const fmpz *a,
                                         uint64_t len,
                                         int normalisation_of_input);
Lerror_t Lfunc_use_dirichlet_coeffs_acb(Lfunc_t L,
                                        acb_srcptr a,
                                        uint64_t len,
                                        int normalisation_of_input);
```

These front-ends supply Dirichlet coefficients directly, indexed
`a[0] = a_1`, ..., `a[len - 1] = a_len`. They overwrite the coefficient array
instead of multiplying in per-prime factors, so they cannot be combined with
any Euler-factor route or called twice; such conflicts are fatal
{c:macro}`ERR_SUPPLY_CONFLICT`.

The `normalisation_of_input` selector is mandatory:

| Selector | Meaning |
| --- | --- |
| {c:macro}`LFUNC_ALGEBRAIC_NORM` | supplied `a_n` are algebraic; the library applies the `n^{-normalisation}` shift before computing |
| {c:macro}`LFUNC_ANALYTIC_NORM` | supplied `a_n` are already analytic; the library applies no shift |

Any other selector is fatal {c:macro}`ERR_BAD_NORM`. The first coefficient must
be supplied and equal to 1: `len == 0` or `a_1 != 1` is fatal
{c:macro}`ERR_A1_NOT_ONE`. A positive-length `NULL` coefficient array is fatal
{c:macro}`ERR_BAD_SUPPLY`.

The `fmpz` form takes exact integer coefficients. The `acb` form trusts the
supplied balls; for `a_1`, the real ball must contain 1 and the imaginary ball
must contain 0. After shifting to analytic normalisation, each coefficient is
checked against the degree's Euler-product bound; a definite violation is
fatal {c:macro}`ERR_COEFF_BOUND`, often indicating the wrong normalisation
selector or corrupted input.

If `len < nmax`, the library reduces `nmax` to `len` and returns the warning
bit {c:macro}`ERR_INSUFF_EULER`; surplus coefficients are ignored. In the
default TURING build, raw coefficients still have enough data for the Turing RH
check. A Buthe-only build reports {c:macro}`ERR_RH_UNAVAILABLE` for raw
coefficient input because Buthe's method needs per-prime data.

### Deliberately reducing the bound: `Lfunc_reduce_nmax`

```c
bool Lfunc_reduce_nmax(Lfunc_t L, uint64_t nmax);
```

If you cannot reach the original {c:func}`Lfunc_nmax`, call this to declare a
smaller trusted bound. For push-style supply, call it after the last successful
push and before {c:func}`Lfunc_compute`. For callback, factor-array, and raw
coefficient routes, call it before the supply front-end if you want deliberate
truncation without {c:macro}`ERR_INSUFF_EULER`; once a short supply call has
already returned that warning, a later reduction cannot erase it. The function
returns `true` only when it actually lowered the current bound (`nmax` strictly
below the current `M`). **The library takes your word for it and does not
check** that the reduced bound still yields a correct result, so use it only
when you understand the consequence.

## Querying the results

Call these only after {c:func}`Lfunc_compute` has returned without a fatal
flag. The `arb_srcptr` / `acb_srcptr` return values are **borrowed**: they point
into the `Lfunc_t` and stay valid until {c:func}`Lfunc_clear`. Do not free
them; copy out (`arb_set` / `acb_set`) anything you need to outlive the object.

### Rank: `Lfunc_rank`

{c:func}`Lfunc_rank` returns the analytic rank (order of vanishing at the
central point) as an `int64_t`. A returned rank of 0 or 1 is rigorous only when
there was no short supply and the RH check ran and succeeded: there must be no
{c:macro}`ERR_INSUFF_EULER`, {c:macro}`ERR_RH_ERROR`,
{c:macro}`ERR_RH_UNAVAILABLE`, {c:macro}`ERR_NO_RANK`, or
{c:macro}`ERR_CONFLICT_RANK` warning. A value greater than 1 is not certified.
If the rank could not be determined you get {c:macro}`ERR_NO_RANK`; if it
disagrees with a rank supplied through {c:struct}`Lparams_t` you get
{c:macro}`ERR_CONFLICT_RANK`.

### Root number: `Lfunc_sign` and `Lfunc_sqrt_sign`

{c:func}`Lfunc_sign` returns the sign (root number) `eps` of the functional
equation, as a borrowed `acb_srcptr`. It is a complex number of absolute value
1; `|eps| = 1` is a useful sanity check on a computation.

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
carries the same rigour caveat as the rank. To evaluate the L-function exactly
at the central point, use this rather than {c:func}`Lfunc_special_value`.

### Zeros: `Lfunc_zeros`

```c
arb_srcptr Lfunc_zeros(Lfunc_t L, uint64_t side);
```

Returns the array of imaginary parts of the zeros on the critical line, as
borrowed `arb_t` balls. `side = 0` gives the zeros of `L`; `side = 1` gives
those of the dual L-function (for a self-dual L-function the two agree). When
the rank is 0 or 1, no short-supply/rank warning was raised, and the RH check
ran and succeeded (no {c:macro}`ERR_INSUFF_EULER`,
{c:macro}`ERR_RH_ERROR`, {c:macro}`ERR_RH_UNAVAILABLE`,
{c:macro}`ERR_NO_RANK`, or {c:macro}`ERR_CONFLICT_RANK` warning), the list is
complete up to the height reached; otherwise some zeros may be missing. At most
{c:macro}`MAX_ZEROS` (256) zeros are stored per side. If the zeros could not be
refined to the target precision you get the {c:macro}`ERR_ZERO_PREC` warning.

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
lets you choose between `L` and the completed `Lambda`.

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
| {c:macro}`ERR_SUPPLY_CONFLICT` | incompatible, duplicate, or out-of-order supply route/factor |
| {c:macro}`ERR_A1_NOT_ONE` | raw Dirichlet-coefficient supply did not provide `a_1 = 1` |
| {c:macro}`ERR_COEFF_BOUND` | a supplied raw coefficient definitely exceeds the degree's Euler-product bound |
| {c:macro}`ERR_BAD_NORM` | invalid `normalisation_of_input`; use {c:macro}`LFUNC_ALGEBRAIC_NORM` or {c:macro}`LFUNC_ANALYTIC_NORM` |
| {c:macro}`ERR_BAD_SUPPLY` | invalid supply argument, including a positive-length `NULL` factor or coefficient array |

### Warning codes

| Macro | Meaning |
| --- | --- |
| {c:macro}`ERR_SOME_DATA` | had sensible data, but not all the way through the Turing zone |
| {c:macro}`ERR_ZERO_PREC` | could not isolate zeros to the target precision |
| {c:macro}`ERR_RH_ERROR` | the RH check on the computed zeros failed |
| {c:macro}`ERR_INSUFF_EULER` | a factor or coefficient supply ended before the expected bound, and `nmax` was reduced |
| {c:macro}`ERR_NO_RANK` | could not determine the rank |
| {c:macro}`ERR_CONFLICT_RANK` | the computed rank disagreed with the supplied one |
| {c:macro}`ERR_DBL_ZERO` | a stationary point failed to converge (a double zero?) |
| {c:macro}`ERR_SPEC_PREC` | could not reach the target error bound in a special value |
| {c:macro}`ERR_G_OUTFILE` | could not open the file to write the G-data cache |
| {c:macro}`ERR_RH_UNAVAILABLE` | RH verification was skipped or unavailable for this build/supply route |

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
so the library caches it to disk. The filename is built from the sorted analytic
shifts (each `mus[i] + normalisation`), one `_<value>` field per shift formatted
to one decimal place: a degree-2 object with analytic shifts `[0.5, 1.5]` is
cached as `g_0.5_1.5`. The cache is only consulted in the default mode (you left
`gprec` at 0, kept the {c:macro}`DEFAULT_TARGET_PREC` target, and set a
`cache_dir`).

Each cache file begins with a self-describing `GCACHE` header recording a format
version, the degree, the `gprec` the grid was computed at, and the sorted shifts.
On a subsequent run the library reads that header and **validates it against the
current request**:

- If the header matches (same version, degree, shifts) and the grid was computed
  at sufficient precision, the body is read back and used as-is.
- If it does not match, or was computed at too low a precision (a *stale* cache),
  the library recomputes the G data and **overwrites** the file rather than
  trusting or aborting on it. A foreign or older-format file at the same name is
  handled the same way.
- If the header is valid but the body cannot be parsed (a *corrupt* file), the
  run fails with the fatal {c:macro}`ERR_G_INFILE`.

{c:member}`cache_dir <Lparams_t.cache_dir>` (settable through
{c:func}`Lfunc_init_advanced`) chooses the directory for these files; left
`NULL` the cache is disabled, and {c:func}`Lfunc_init` defaults it to the current
working directory.

```{note}
Because the filename encodes only the shifts (not the degree, precision, or
version), two requests with the same analytic shifts map to the same filename;
the `GCACHE` header is what actually distinguishes them and triggers a recompute
on a mismatch. For a hermetic computation, still point {c:member}`cache_dir
<Lparams_t.cache_dir>` at a clean directory or remove any `g_*` files first: the
header validation prevents a *stale* cache from being trusted, but a clean
directory avoids the recompute-and-overwrite churn entirely. A grid that does
not extend far enough for the conductor raises {c:macro}`ERR_G_EXTENT`.
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
  {c:func}`Lfunc_use_lpolys_fmpz` consecutive-prime factor-array route.
- [`tau_dirichlet.cpp`](https://github.com/edgarcosta/lfunctions/blob/main/examples/tau_dirichlet.cpp):
  the Ramanujan tau example supplied by raw coefficients through
  {c:func}`Lfunc_use_dirichlet_coeffs_fmpz`, demonstrating
  {c:macro}`LFUNC_ALGEBRAIC_NORM` and the default Turing RH check on raw
  coefficient input.

## Full header reference

The declarations below are generated from
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h)
by Doxygen and rendered through Breathe, so the header stays the single source
of truth for the per-function reference.

```{doxygenfile} glfunc.h
:project: liblfun
```
