# LMFDB regression curves for the smalljac -> rational.exe pipeline

A small, hand-curated set of curves drawn from the
[LMFDB](https://www.lmfdb.org), with certified reference values, intended to
drive an end-to-end CI regression of `examples/rational.c` (`rational.exe`)
once the smalljac point-counting step is wired into CI (bead `lfunctions-4jd`).

Files in this directory:

| file               | role                                                                 |
|--------------------|----------------------------------------------------------------------|
| `curves.smalljac`  | the curve list, in **smalljac's input format** (one curve per line)  |
| `reference.txt`    | parallel **certified reference values** (rank, root number, leading Taylor coeff, first zero) |
| `README.md`        | this file: provenance, selection rationale, reproduction, tolerances |

This is **data only**: nothing here is built or linked. No source under `src/`
or `examples/` is modified by this bead.

## The pipeline this feeds

```
curves.smalljac                    (this dir, smalljac INPUT format)
      |  smalljac driver: count points -> good Euler factors,
      |  splice in the hard-coded bad factors from each line
      v
euler-factor file                  (rational.exe INPUT format:
      |                             label:degree:conductor:weight:[mus]:[[euler_factors]])
      v
rational.exe  --->  Rank, Epsilon, leading Taylor coeff, first zero
      |
      v
compare against reference.txt      (the CI assertion; bead lfunctions-4jd)
```

`rational.exe`'s two input formats are documented at the top of
`examples/rational.c`. The **smalljac** input format used by `curves.smalljac`
is documented at the top of the (since-removed) `examples/smalljac.cpp`; the
authoritative field-parsing loop is in that file at git commit `660d3c8~1`.

### smalljac input line format

```
label:symdegree:conductor:curve_str:bad_euler_factors
```

- `label` — LMFDB label.
- `symdegree` — symmetric-power degree; **always `1`** here (we want `L(C)`,
  not `Sym^k L(C)`). smalljac uses this field to select the symmetric-power
  Euler-factor formulas (only implemented for `E/Q`).
- `conductor` — the arithmetic conductor of the L-function.
- `curve_str` — the defining equation:
  - elliptic curve: `[a1,a2,a3,a4,a6]` (the five a-invariants);
  - genus 2 curve: `[f(x), h(x)]` for the model `y^2 + h(x) y = f(x)`.
- `bad_euler_factors` — `[[p1,p2,...],[Lp1,Lp2,...]]`: the bad primes and their
  local L-factors `Lp(T) = sum_i coeff_i T^i` (coefficients low -> high degree).
  smalljac computes the *good* factors itself; the bad factors must be supplied.

The resulting `mus` and `weight` are not on the line: smalljac derives
`degree = 2*genus`, `mus = [0,...,0, 1,...,1]` (half zeros, half ones), and
`weight = symdegree*0.5`. This matches the `rational.c` examples
`11a2:2:11:1:[0,1]:...` (EC) and `169.a.169.1:4:169:1:[0,0,1,1]:...` (genus 2).

## Selection rationale / coverage

Chosen to satisfy bead `lfunctions-6dv`'s coverage targets while keeping every
object cheap enough for CI (smallest available conductor in each bucket; the
deg-2 set reuses the canonical objects already asserted by the `ec_*.cpp`
examples, so the smalljac path is cross-checked against the hand-written path).

**Elliptic curves over Q — distinct analytic ranks 0..4** (degree-2 L-functions):

| label        | conductor | analytic rank | notes                                   |
|--------------|-----------|---------------|-----------------------------------------|
| `11.a2`      | 11        | 0             | smallest conductor; matches `examples/`-style asserts |
| `37.a1`      | 37        | 1             | also `examples/ec_37.a1.cpp`            |
| `389.a1`     | 389       | 2             | also `examples/ec_389.a1.cpp`           |
| `5077.a1`    | 5077      | 3             | also `examples/ec_5077.a1.cpp`; named in the bead |
| `234446.a1`  | 234446    | 4             | rank-4 example named in the bead        |

**Genus 2 curves over Q — distinct analytic ranks 0..3** (degree-4 L-functions),
all geometrically simple with `end_alg = Q` and Sato-Tate `USp(4)` (i.e.
primitive degree-4 L-functions, the generic case):

| label             | conductor | analytic rank | rank proved? |
|-------------------|-----------|---------------|--------------|
| `249.a.249.1`     | 249       | 0             | no (upper bound, tight) |
| `587.a.587.1`     | 587       | 1             | yes          |
| `3319.a.3319.1`   | 3319      | 2             | no (upper bound, tight) |
| `25913.a.25913.1` | 25913     | 3             | no (upper bound, tight) |

`587.a.587.1` is the smallest-conductor analytic-rank-1 genus 2 curve in LMFDB;
`3319.a.3319.1` and `25913.a.25913.1` are the smallest at ranks 2 and 3. The
genus 2 analytic ranks above 1 are upper bounds in LMFDB
(`analytic_rank_proved = false`), which is expected: it does not affect the
numeric L-data the library computes, only the BSD-side interpretation.

**Square case `L(C) = L(E)^2`** (degree-4, but reducible):

| label             | conductor | rank | end_alg   | st_group | why                                  |
|-------------------|-----------|------|-----------|----------|--------------------------------------|
| `196.a.21952.1`   | 196       | 0    | `M_2(Q)`  | `E_1`    | Jac ~ E x E with E,E isogenous; the L-function L-label is `4-14e2-1.1-c1e2-0-0` and **every zero has multiplicity 2** |

This is the confirmed example from the bead (`end_alg = M_2(Q)`,
`is_simple_geom = false`). The doubled zeros make it a strong structural test
of the zero-finding / rank logic. More such cases can be found with the filter
`end_alg LIKE 'M_2%' AND is_simple_geom = false`.

**Product case `L(C) = L(E1) L(E2)`, E1 not isogenous to E2** (the optional
stretch target):

| label          | conductor | rank | end_alg   | st_group        | why                                   |
|----------------|-----------|------|-----------|-----------------|---------------------------------------|
| `294.a.294.1`  | 294       | 0    | `Q x Q`   | `SU(2)xSU(2)`   | Jac splits as a product of two non-isogenous elliptic curves; zeros interleave two distinct EC L-functions |

Found with `end_alg LIKE 'Q x Q'` ordered by conductor.

## Reproducibility — exact LMFDB queries

The mirror is read-only PostgreSQL; table/column names are case-sensitive
(`"Lhash"` must be quoted). All values below were pulled 2026-06-02.

### 1. Picking the curves

Elliptic curves (canonical small-conductor representative of each rank):

```sql
-- rank r representative used here; r in 0..4
SELECT lmfdb_label, conductor, analytic_rank, ainvs, bad_primes
FROM ec_curvedata
WHERE lmfdb_label IN ('11.a2','37.a1','389.a1','5077.a1','234446.a1')
ORDER BY conductor;
```

Genus 2, generic, smallest conductor per rank:

```sql
SELECT DISTINCT ON (analytic_rank)
       analytic_rank, label, cond, analytic_rank_proved, root_number,
       end_alg, is_simple_geom, st_group, eqn, bad_primes, bad_lfactors
FROM g2c_curves
WHERE analytic_rank IN (0,1,2,3) AND end_alg = 'Q'
ORDER BY analytic_rank, cond ASC;
```

Square case and product case:

```sql
SELECT label, cond, analytic_rank, end_alg, is_simple_geom, st_group, eqn,
       bad_primes, bad_lfactors
FROM g2c_curves
WHERE label = '196.a.21952.1';                 -- square, L = L(E)^2

SELECT label, cond, analytic_rank, end_alg, is_simple_geom, st_group, eqn,
       bad_primes, bad_lfactors
FROM g2c_curves
WHERE end_alg LIKE 'Q x Q'
ORDER BY cond ASC LIMIT 5;                       -- product; we took 294.a.294.1
```

### 2. The reference values

Elliptic-curve L-data (join via the object URL):

```sql
SELECT i.url, i.label AS lfunc_label,
       l.order_of_vanishing, l.root_number, l.leading_term,
       l.accuracy, l.gamma_factors, l.positive_zeros
FROM lfunc_instances i
JOIN lfunc_lfunctions l ON i."Lhash" = l."Lhash"
WHERE i.url IN ('EllipticCurve/Q/11/a','EllipticCurve/Q/37/a',
                'EllipticCurve/Q/389/a','EllipticCurve/Q/5077/a',
                'EllipticCurve/Q/234446/a')
ORDER BY l.conductor;
```

Genus 2 L-data (join via the genus-2 curve's own `Lhash`):

```sql
SELECT g.label, l.label AS lfunc_label,
       l.degree, l.order_of_vanishing, l.root_number, l.leading_term,
       l.accuracy, l.gamma_factors, l.positive_zeros
FROM g2c_curves g
JOIN lfunc_lfunctions l ON g."Lhash" = l."Lhash"
WHERE g.label IN ('249.a.249.1','587.a.587.1','3319.a.3319.1',
                  '25913.a.25913.1','196.a.21952.1','294.a.294.1')
ORDER BY g.cond;
```

`reference.txt` records `order_of_vanishing` -> rank, `root_number`,
`leading_term`, and `positive_zeros[0]` -> first_zero from these queries.

### 3. Deriving the smalljac line from the LMFDB row

- **EC `curve_str`** = the `ainvs` array `[a1,a2,a3,a4,a6]` verbatim.
- **EC bad factor** at the (single, here) bad prime p, from
  `ec_localdata.reduction_type`:
  `+1` split -> `[1,-1]`; `-1` non-split -> `[1,1]`; `0` additive -> `[1]`.
  (`SELECT lmfdb_label, prime, reduction_type FROM ec_localdata WHERE ...`)
- **Genus 2 `curve_str`** = `[f(x), h(x)]` built from
  `g2c_curves.eqn = [[f_coeffs],[h_coeffs]]` (coeffs low -> high degree),
  e.g. `[[0,-2,-1,1],[1,1,0,1]]` -> `[-2*x - x^2 + x^3, 1 + x + x^3]`.
  Term order / spacing is irrelevant — smalljac parses the polynomial string.
- **Genus 2 bad factors** = transpose of `g2c_curves.bad_lfactors`, which is
  `[[p,[Lp_coeffs]],...]`, into smalljac's `[[p,...],[Lp,...]]`.

## Tolerances

The library returns **certified `arb`/`acb` balls**; the strings here are the
*reference* to test those balls against, so the relevant question is how many
digits of each reference string are trustworthy.

- **Elliptic curves**: LMFDB `accuracy = 100` bits for these L-functions
  (`SELECT label, accuracy FROM lfunc_lfunctions WHERE label IN (...)`). The
  zero/value strings stored carry ~16 significant figures, all correct. Assert
  with `arb_contains` (the computed ball should contain the reference) or check
  `|computed - reference| < 1e-12`.
- **Genus 2 curves**: `accuracy` is NULL; `positive_zeros`/`leading_term` are
  stored as ~15-16 digit floats. Use `1e-11`.

For the **root number**, prefer asserting the mathematical invariant
`|Lfunc_sign(L)| = 1` (always true, LMFDB-independent) in addition to matching
the recorded sign.

## TODO / open items for lfunctions-4jd

All reference values above are present (none left as TODO). Items that belong to
the *consuming* CI bead, not this data bead:

- Wire `curves.smalljac` through the smalljac driver to emit the
  `rational.exe` input file, then assert `rational.exe` output against
  `reference.txt`. `rational.c` currently only *prints* results (its `main`
  has a `// TODO write output to output`), so lfunctions-4jd will either parse
  stdout or extend `rational.c` to emit a comparable record.
- Optionally add a non-central special value `L(s0)` per curve to
  `reference.txt` (the `ec_*.cpp` examples assert `L(1.5)` and `L(2.5)` in
  analytic normalization). These are not stored verbatim in
  `lfunc_lfunctions`, so they would be generated/recorded at wiring time;
  hence omitted here to keep every value LMFDB-sourced and reproducible.
