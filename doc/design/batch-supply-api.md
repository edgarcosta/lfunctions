# Batch supply API for Dirichlet coefficients and Euler factors

Status: **design / decision record** (bead `lfunctions-wor`). No code is changed by
this document. It specifies the API to be implemented by a follow-up bead and
records the resolution of the seven open questions raised in the bead.

## 1. Motivation and current state

Today the only way to feed arithmetic data into an `Lfunc_t` is through Euler
factors, in one of two forms (`include/glfunc.h`):

- `Lfunc_use_all_lpolys(L, callback, param)` — the library drives a prime sieve
  over `p <= Lfunc_nmax(L)` and calls back for each prime; returning a zero
  polynomial short-circuits and reduces `nmax`.
- `Lfunc_use_lpoly(L, p, poly)` — the caller pushes one `acb_poly_t` factor.

Both funnel through the internal `use_lpoly()` (`src/coeff.c:121`), which

1. normalises the factor by multiplying the coefficient of `T^m` by
   `p^{-m·normalisation}` (`src/coeff.c:135-145`) to move from the algebraic to
   the analytic normalisation;
2. inverts the local series to the needed length with `acb_poly_inv_series`
   (`src/coeff.c:150`);
3. hands the inverse (and, under `BUTHE`, the forward factor) to
   `use_inv_lpoly()` (`src/coeff.c:88-119`), which multiplicatively sieves the
   inverse-series coefficients into the global Dirichlet array `L->ans`
   (initialised to all-ones in `Lfunc_nmax`, `src/coeff.c:75-76`) and, under
   `BUTHE`, accumulates the Buthe `Wf` sum via `wf()`.

Two gaps motivate this bead:

- **No way to supply `a_n` directly.** Callers that already hold Dirichlet
  coefficients (e.g. from a database row, or a non-Euler-product construction)
  must re-package them as local factors.
- **`acb_poly_t` is the only accepted factor type**, and the only batch entry
  point is a C function pointer. Function pointers bind poorly to Magma and most
  other host languages (the planned port), and integer data (exact `fmpz` /
  `fmpz_poly`) must be hand-converted to `acb` before it can be passed.

Goal: add ergonomic, data-driven, binding-friendly **array** front-ends for the
common case, while **retaining** `Lfunc_use_all_lpolys` as the streaming/lazy
backend for large `nmax`, and routing every new entry point through the *same*
internal fold so there is no second copy of the coefficient math or the
short-circuit logic.

## 2. Proposed public surface

Four new entry points, plus one shared internal length-control helper. Names are
final recommendations; the bead listed them as illustrative.

```c
// --- Raw Dirichlet coefficients a_n, indexed n = 1..len (a[0] = a_1) -------
// `normalisation_of_input` says which normalisation the supplied a_n are in:
//   ALGEBRAIC_NORM  -> library applies the n^{-normalisation} shift,
//   ANALYTIC_NORM   -> a_n are already analytic; no shift applied.
Lerror_t Lfunc_use_dirichlet_coeffs_fmpz(Lfunc_t L,
                                         const fmpz *a, uint64_t len,
                                         int normalisation_of_input);

Lerror_t Lfunc_use_dirichlet_coeffs_acb(Lfunc_t L,
                                        acb_srcptr a, uint64_t len,
                                        int normalisation_of_input);

// --- Euler factors, one per prime p_k <= nmax, in increasing-prime order ---
// f[k] is the local factor at the k-th prime (f[0] at p=2). `len` is the count
// of factors supplied; the library sieves primes internally and consumes the
// k-th factor for the k-th prime. Factors are in algebraic normalisation
// (same convention as Lfunc_use_lpoly); the p^{-m·normalisation} shift is
// applied by the shared use_lpoly path.
Lerror_t Lfunc_use_lpolys_acb(Lfunc_t L,
                              const acb_poly_struct *f, uint64_t len);

Lerror_t Lfunc_use_lpolys_fmpz(Lfunc_t L,
                               const fmpz_poly_struct *f, uint64_t len);
```

Supporting constants (new, in `glfunc.h`), used for the
`normalisation_of_input` argument so the contract is explicit rather than a bare
`int` or a guessed default:

```c
#define ALGEBRAIC_NORM (0) // supplied a_n are algebraic; apply n^{-normalisation}
#define ANALYTIC_NORM  (1) // supplied a_n already analytic; no shift
```

Rationale for the per-input normalisation flag rather than always-algebraic:
unlike Euler factors (which are essentially always handed over in the algebraic
normalisation, matching `Lfunc_use_lpoly`), raw `a_n` are frequently already in
the analytic normalisation when they come out of another tool. Forcing one
convention would make the other caller silently wrong. The Euler-factor
front-ends keep the single algebraic convention to stay bug-for-bug compatible
with `Lfunc_use_lpoly`; we do **not** add a flag there.

### 2.1 `acb_poly_struct *` vs `acb_poly_t *`

`acb_poly_t` is `acb_poly_struct[1]`, so an array of polynomials is naturally an
`acb_poly_struct *` (a contiguous block of structs), which is exactly what
`f[k]` indexes and what FLINT's own vector code uses. The signatures take
`const acb_poly_struct *` / `const fmpz_poly_struct *` for that reason. A host
binding constructs the block once and passes the base pointer — no array of
pointers, no per-element marshalling.

## 3. Funnelling into the shared internal path (Q6, Q2)

The single most important constraint: **all five supply routes (callback, push,
two coefficient arrays, two factor arrays) end at the same internal fold.** The
implementation refactors `src/coeff.c` so that:

```
                         Lfunc_use_all_lpolys (callback, streaming)  ─┐
                         Lfunc_use_lpoly (push acb_poly)             ─┤
   Lfunc_use_lpolys_acb  (array of acb_poly)  ── per p ───────────────┤──> use_lpoly(L, p, acb_poly)
   Lfunc_use_lpolys_fmpz (array of fmpz_poly) ── convert+per p ───────┘        │
                                                                               │ normalise (p^{-m·norm}),
                                                                               │ inv_series, then
                                                                               ▼
                                                                       use_inv_lpoly(L, p, ...) ─> multiplicative sieve into L->ans
   Lfunc_use_dirichlet_coeffs_fmpz ─┐
   Lfunc_use_dirichlet_coeffs_acb  ─┴── normalise (n^{-norm}) ─────────────────────────────────> write directly into L->ans
```

Concretely:

- **Euler-factor arrays** are *thin*. `Lfunc_use_lpolys_acb` iterates a
  primesieve exactly like `Lfunc_use_all_lpolys`, and for prime `p_k` calls the
  existing `use_lpoly(L, p_k, f[k])`. Nothing about normalisation, inversion, or
  the sieve is duplicated. `Lfunc_use_lpolys_fmpz` is identical except it first
  converts `f[k]` (an `fmpz_poly`) into a scratch `acb_poly` via
  `acb_poly_set_fmpz_poly` (exact, no rounding) and passes that to `use_lpoly`.
  Buthe/`wf` therefore keeps working unchanged for both, because real per-prime
  forward *and* inverse factors still exist.

- **Raw `a_n` arrays** cannot go through `use_lpoly` (there is no per-prime
  factor to invert). They write straight into `L->ans[n-1]`. To avoid a second
  copy of *anything*, the implementation must:

  - reuse `Lfunc_nmax` to size/allocate/initialise `L->ans` (it already runs the
    `L->M` derivation and the `acb_set_ui(L->ans[i],1)` initialisation; the
    raw-coefficient path **overwrites** `ans` rather than multiplying into it, so
    it must run *instead of*, not after, any Euler supply — see §6 and OQ-3);
  - factor the **insufficient-supply** decision into one shared helper, next.

### 3.1 One shared short-circuit helper (Q2)

The callback path currently inlines the "ran out of factors" logic
(`src/coeff.c:198-206`): on a zero poly it shrinks `L->M = p-1`, shrinks
`buthe_M` under `BUTHE`, and OR-s in `ERR_INSUFF_EULER`. The push path
(`Lfunc_use_lpoly`) has no such guard, and `Lfunc_reduce_nmax` duplicates a
similar `M`/`buthe_M` clamp (`src/coeff.c:219-232`).

Design decision: extract a single internal helper

```c
// reduce the working coefficient count to `new_M` (because we ran out of data
// at, or were told to stop at, index/prime new_M+1), clamp buthe_M, and return
// the warning bit to OR into the caller's accumulator.
static Lerror_t shrink_M(Lfunc *L, uint64_t new_M, bool insufficient);
```

and route the callback short-circuit, `Lfunc_reduce_nmax`, and the array
length-shortfall path through it. `insufficient == true` sets
`ERR_INSUFF_EULER`; `Lfunc_reduce_nmax` passes `false` (an explicit, trusted
reduction is not an error). This removes the existing `M`/`buthe_M` clamp
duplication as a side benefit.

## 4. Indexing and length (Q1)

Resolved as follows.

- **Raw `a_n`**: indexed `n = 1..len`, so `a[0] == a_1`. The library reads
  `min(len, nmax)` of them where `nmax = Lfunc_nmax(L)`. There is no
  "0th coefficient." If `len < nmax`, see §5.

  - Convention check: `a_1` must be `1`. If `a[0] != 1` (exactly, for `fmpz`; or
    the `acb` ball failing to contain 1, see OQ-5), OR in a fatal
    `ERR_A1_NOT_ONE` (see §7/§10). This catches off-by-one packing (caller
    passed `a_0..a_{len-1}`) immediately.

- **Euler-factor arrays**: indexed over **consecutive primes**, `f[k]` for the
  k-th prime (`f[0]` at `p = 2`), length `len = pi(nmax)` for a complete supply.
  The library generates the primes itself (it already owns a primesieve);
  **the caller does not pass the primes.**

  Rationale for *not* also taking a `primes[]` array: the set of primes consumed
  is fully determined by `nmax`, which the library computes. Passing primes
  invites a mismatch (caller's prime list disagreeing with the internal sieve)
  with no upside, and doubles the marshalling burden for bindings. A caller who
  genuinely wants to control which primes are used should use the push API
  (`Lfunc_use_lpoly`) or the callback. (Open question OQ-1 in §9 revisits a
  sparse/`(p, poly)`-pairs variant; the recommendation is to defer it.)

## 5. Insufficient-supply semantics (Q2)

If a supplied array is shorter than the library needs, the behaviour mirrors the
callback's existing zero-poly short-circuit, via the shared `shrink_M` helper
(§3.1):

- **Euler-factor array** with `len < pi(nmax)`: factors are consumed for the
  first `len` primes; on running out at prime `p`, `shrink_M(L, p-1, true)`
  reduces `M` (and `buthe_M`) and sets `ERR_INSUFF_EULER`. Identical outcome to
  the callback returning a zero poly at `p`. This is a **warning**, not fatal:
  the computation proceeds with the reduced `M`, exactly as today.

- **Raw `a_n` array** with `len < nmax`: coefficients are written for
  `n = 1..len`; then `shrink_M(L, len, true)` reduces `M` to `len` and sets
  `ERR_INSUFF_EULER`. (Index `len` is a coefficient index here, not a prime,
  but the clamp semantics are the same.)

- Supplying **more** than needed (`len > nmax` / `len > pi(nmax)`) is allowed
  and not an error: the surplus is ignored. This matches the spirit of
  `Lfunc_reduce_nmax` ("takes your word for it") and lets a caller hand over a
  fixed-size buffer.

`ERR_INSUFF_EULER` already exists and already has an `fprint_errors` branch
(`src/glfunc.c:38`); reuse it for the array shortfall. Its message text may want
broadening from "Euler" to "Euler factors / coefficients" so it reads correctly
for the raw-`a_n` case (minor, noted for the implementation bead).

## 6. Normalisation contract (Q3)

The library runs entirely in the **analytic** normalisation
(`L->mus[i] = mus[i] + normalisation`, coefficients centred at `s = 1/2`). The
shift between algebraic and analytic is `n^{-normalisation}` on the coefficient
`a_n` — equivalently `p^{-m·normalisation}` on the `T^m` coefficient of the
local factor at `p`, which is exactly what `use_lpoly` applies
(`src/coeff.c:135-145`).

Per input form:

- **Euler factors** (`acb_poly` and `fmpz_poly` arrays): same contract as the
  existing `Lfunc_use_lpoly` — factors are in the **algebraic** normalisation,
  and the library applies `p^{-m·normalisation}`. The `fmpz_poly` form converts
  exactly to `acb_poly` first (§3), then takes the identical shift. No flag.

- **Raw `a_n`** (`fmpz` and `acb` arrays): the caller declares the input
  normalisation via `normalisation_of_input`.
  - `ALGEBRAIC_NORM`: the library multiplies `a_n` by `n^{-normalisation}`
    (the direct analogue of the per-factor shift; computed as
    `exp(-normalisation·log n)` at working precision, matching the `arb_exp` of
    `-normalisation·log p` in `use_lpoly`).
  - `ANALYTIC_NORM`: no shift; the `a_n` are written as given.
  - When `normalisation == 0` the two are identical, but the flag is still
    required (no silent default) so the intent is explicit and auditable.

Documentation must state this contract on each function and cross-reference the
existing `mus`/`normalisation` discussion in `glfunc.h`.

## 7. Buthe / RH behaviour (Q4)

Buthe's method (compiled in only under `#define BUTHE`; `TURING` is the default
RH check) needs, per prime, *both* the forward local factor `fp` and its inverse
`fp1` — `wf()` (`src/buthe.c:139`) multiplies coefficients of the two together.

- **Euler-factor arrays** (`acb_poly`/`fmpz_poly`): per-prime forward and
  inverse factors exist (the `fmpz_poly` form is converted to `acb_poly` first),
  so `wf()` runs exactly as for `Lfunc_use_lpoly` / the callback. **No change to
  RH behaviour.**

- **Raw `a_n`**: there are **no per-prime factors**, so `wf()` cannot be fed and
  the Buthe sum cannot be formed. The library must **disable RH verification and
  say so** — never silently skip it. Introduce a new **warning** bit (upper 32
  of `Lerror_t`):

  ```c
  #define ERR_RH_UNAVAILABLE ((uint64_t) 1<<41) // RH check skipped: no per-prime
                                                 // factors (raw a_n supplied)
  ```

  (`1<<40` is the last warning currently in use, `ERR_G_OUTFILE`; `1<<41` is the
  next free upper bit.) Set this bit when a raw-`a_n` supply path is used and an
  RH check would otherwise have run, and add a matching `fprint_errors` branch in
  `src/glfunc.c`:

  > "RH check skipped: raw Dirichlet coefficients were supplied, so no
  > per-prime Euler factors are available for the Buthe/Turing verification."

  Mechanically: the raw-`a_n` front-ends set a flag on the `Lfunc` struct (e.g.
  `L->no_lpolys = true`); `Lfunc_compute` consults it and, instead of calling
  the RH check, OR-s `ERR_RH_UNAVAILABLE` into its return. Under `BUTHE`,
  `buthe_M`/`buthe_Wf` are simply never populated; under `TURING`, the
  Turing check is likewise skipped. The downstream effect is the same as today
  when an RH check cannot complete: `Lfunc_zeros` may be incomplete and the
  caller is warned via the error bitfield. (Note: per the project's known issue,
  the Turing RH count is currently miscalibrated for degree >= 3 and already
  raises `ERR_RH_ERROR`; `ERR_RH_UNAVAILABLE` is a distinct, intentional "we
  never even tried" signal for the raw-`a_n` case and must not be conflated with
  it.)

- **Out of scope (recommended):** reformulating an RH check directly from `a_n`
  (without per-prime factors) is a research task, not an API task. The bead
  flags it as "likely out"; this design **defers** it (OQ-2 in §9). Until then,
  callers who need a verified zero set must supply Euler factors.

## 8. Type design and certified discipline (Q5)

Consistent with the project's certified-ball direction (every numeric input is a
ball or an exact integer; no plain `double` coefficient inputs):

- **`acb` / `acb_poly` inputs are balls and are trusted as given.** The library
  does not re-verify them; the caller asserts the ball encloses the true value.
  This matches how `Lfunc_use_lpoly` already treats its `acb_poly_t`.

- **`fmpz` / `fmpz_poly` inputs are exact integers** and are converted to
  **exact** `acb` / `acb_poly` (zero-radius balls) via `acb_set_fmpz` /
  `acb_poly_set_fmpz_poly`. This is the *preferred* form whenever the data is
  genuinely integral (most arithmetic L-functions of interest), because it
  carries no rounding and binds trivially from languages with big integers
  (Magma `RngIntElt`, FLINT `fmpz`).

- **No plain-float entry points.** A caller holding `double`s must build an
  `acb` (and own the rounding decision) before calling. This is deliberate and
  preserves the end-to-end certified interface.

## 9. Open questions and recommendations

- **OQ-1 (sparse / explicit-prime factor supply).** Should there be a variant of
  the Euler-factor array that also takes the primes (e.g. parallel
  `uint64_t primes[]` + `acb_poly_struct factors[]`, or an array of `(p, poly)`
  pairs), for callers who skip primes or supply out of order?
  **Recommendation: defer.** The push API already covers arbitrary `(p, poly)`,
  and the consecutive-prime array covers the common dense case. Adding a sparse
  array now multiplies the surface for a thin benefit. Revisit if a concrete
  binding needs it.

- **OQ-2 (RH from raw `a_n`).** Is a Buthe/Turing reformulation that does not
  need per-prime factors feasible? **Recommendation: out of scope here**; track
  as a separate research bead. For now raw-`a_n` supply yields
  `ERR_RH_UNAVAILABLE`.

- **OQ-3 (mixing supply forms).** May a caller mix, e.g., push some factors then
  call an array, or supply raw `a_n` *and* factors? **Recommendation: forbid
  mixing raw-`a_n` with any Euler-factor route, and forbid calling a coefficient
  array more than once.** The raw-`a_n` path *overwrites* `L->ans`, whereas the
  factor paths *multiply into* it (starting from the all-ones init); interleaving
  the two is meaningless. The implementation should track that a supply step has
  occurred and OR a fatal `ERR_SUPPLY_CONFLICT` (new, low-32 bit) if an
  incompatible second supply call is made. Mixing *callback + push +
  factor-array* (all multiplicative) is internally consistent and may be
  permitted, but the recommendation is to keep it simple: **one supply call per
  `Lfunc_t`** in the first implementation, relaxing later only if a use case
  appears.

- **OQ-4 (`a_1 != 1`).** Strictly, a normalised L-function has `a_1 = 1`. The
  design treats `a_1 != 1` as fatal (`ERR_A1_NOT_ONE`, §4/§10). Confirm no
  in-scope object legitimately violates this; if one does, downgrade to a
  warning. (For the factor paths `a_1 = 1` is automatic from the all-ones init
  and the local factors' constant term, so this only concerns the raw-`a_n`
  paths.)

- **OQ-5 (acb `a_1` exactness).** For `Lfunc_use_dirichlet_coeffs_acb`, do we
  require `a[0]` to *contain* 1 (ball overlaps 1) or to be *exactly* 1?
  **Recommendation: require the ball to contain 1** (use `arb_contains_si` on the
  real part and `arb_contains_zero` on the imaginary part); rejecting a
  legitimately wide-but-correct ball would be wrong, and the certified contract
  is "the ball encloses the truth."

## 10. New error/flag codes introduced (summary)

| Code | Bit | Kind | Meaning |
|------|-----|------|---------|
| `ERR_SUPPLY_CONFLICT` | `1<<13` (fatal) | fatal | incompatible/duplicate supply call (see OQ-3) |
| `ERR_A1_NOT_ONE` | `1<<14` (fatal) | fatal | supplied `a_1` is not 1 (raw-`a_n` paths) |
| `ERR_RH_UNAVAILABLE` | `1<<41` (warning) | warning | RH check skipped: raw `a_n`, no per-prime factors |
| `ALGEBRAIC_NORM` / `ANALYTIC_NORM` | n/a | enum-ish `int` | `normalisation_of_input` selector for raw-`a_n` |

Both new error codes need (a) a `#define` in `glfunc.h` in the correct
fatal/warning half, and (b) a `fprint_errors` branch in `src/glfunc.c`
(per the project convention that every flag prints one human line). The next
free *fatal* (lower-32) bits after `ERR_G_EXTENT` (`= 4096 = 1<<12`) are `1<<13`
and `1<<14`; the next free *warning* (upper-32) bit after `ERR_G_OUTFILE`
(`1<<40`) is `1<<41`.

## 11. Follow-up implementation beads and tests

To be filed after this design is reviewed:

1. **Refactor `coeff.c` to extract `shrink_M`** and route the callback
   short-circuit + `Lfunc_reduce_nmax` through it (pure refactor; the existing
   regression suite is the oracle — `make check` must stay green with no
   behaviour change).
2. **Implement the four array front-ends** + the new error/flag codes +
   `fprint_errors` branches.
3. **Tests** (through the public API only, `arb_overlaps`/`acb_overlaps`/
   `arb_contains` asserts — never printed output):
   - **Equivalence across forms.** A canonical object (e.g. the EC `11.a`
     already used in examples) supplied four ways — callback, `acb_poly` array,
     `fmpz_poly` array, and raw `a_n` (`fmpz`) — must produce zeros / `epsilon` /
     `Lfunc_Taylor` that pairwise overlap. The `acb`-coeff form vs `fmpz`-coeff
     form (exact same numbers) must agree to full precision.
   - **Normalisation flag.** Equivalent inputs (`ALGEBRAIC_NORM` with non-zero
     `normalisation` vs the same numbers pre-shifted, supplied `ANALYTIC_NORM`)
     must agree — the analogue of the existing `[0,1],0.5` == `[0.5,1.5],0`
     invariant.
   - **Insufficient supply.** A short array must set `ERR_INSUFF_EULER`, reduce
     `M`, and otherwise behave like the callback's zero-poly short-circuit
     (compare against a callback run truncated at the same point).
   - **RH-unavailable.** A raw-`a_n` run must set `ERR_RH_UNAVAILABLE` (and must
     *not* set it on a factor-array run); `fatal_error` must remain false for it.
   - **Edge guards.** `a_1 != 1` sets `ERR_A1_NOT_ONE` (fatal); a disallowed
     second supply call sets `ERR_SUPPLY_CONFLICT` (fatal).
   - Use a clean `cache_dir` and `rm -f g_*` for hermetic runs (project
     convention; stale `g_<norm>` files poison results). degree >= 3 tests must
     tolerate the known Turing `ERR_RH_ERROR` warning.
