# smalljac `lpdata` output format (for `test/highdeg/gen.py`)

This documents the exact stdout / output-file format of the smalljac `lpdata`
utility, so the high-degree generator (`gen.py`) can parse it. It was determined
by reading the smalljac source (`lpdata.c`, the `dump_lpoly` callback) **and**
empirically by building smalljac v4.1.3 for genus 2 and running it (see below).

## Which smalljac we build, and why

- `test/highdeg/install_smalljac.sh` builds **smalljac v4.1.3** against
  **ff_poly v1.2.7**, configured for **`SMALLJAC_GENUS = 2`**.
- v4.1.3 with `SMALLJAC_GENUS=2` handles **genus 1 AND genus 2 simultaneously**
  (its `smalljac.h` says: *"genus 1 and 2 ... work when SMALLJAC_GENUS is set to
  2, do not change this"*). That is exactly what the high-degree suite needs:
  genus 1 for the elliptic curves and symmetric powers, genus 2 for the
  degree-4 curves. We do **not** need genus 3.
- We deliberately do **not** use smalljac v5.0: it adds a third dependency,
  `hwlpoly` (the Harvey–Sutherland average-polynomial-time library, pulled in via
  `-lhwlpoly`), which is **not publicly downloadable** from the smalljac site,
  and it is only needed for genus 3. v4.1.3's only deps are **ff_poly + GMP**.
- The smalljac/ff_poly sources are **not on GitHub**; they are tarballs on
  Andrew Sutherland's MIT page (`https://math.mit.edu/~drew/`). The install
  script downloads `smalljac_v4.1.3.tar` and `ff_poly_v1.2.7.tar` from there.

## CLI

```
lpdata file-prefix curve-spec interval/bound [flags jobs job-id]
```

- Output is written to a file named **`<file-prefix>_lpdata.txt`** in the CWD
  (NOT to stdout — stdout only gets progress/summary lines like
  `trace sum is ...`, `Processed N primes ...`, `Output written to file ...`).
- `interval/bound`: a bound `N` means primes in `[1, N]`; `[a..b]` means `[a,b]`.
  Accepts exponential notation, e.g. `10e6`, `2e20`.
- **flags** (4th arg, an integer OR-ed from):
  - `1` = `SMALLJAC_GOOD_ONLY`  — emit only primes of good reduction (what `gen.py` uses).
  - `2` = `SMALLJAC_A1_ONLY`    — emit only the linear coefficient a1.
  - `4` = `SMALLJAC_DEGREE1`    — degree-1 primes only (number-field curves).
  - In v4.1.3 the valid range is `0..7`. (v5.0 adds `8`, `16`, `32`, `64`.)
- `jobs job-id`: parallel sharding; **unused by us** (leave them off → serial).
  If `jobs` is given the output file is `<prefix>_lpdata_<jobs>_<jobid>.txt`.

## Output-file format

First line is a **header** echoing the curve and interval:
```
<curve-spec-in-brackets> <minp> <maxp>
```
e.g. `[0,-1,1,-10,-20] 1 30` or `[x^6+4*x^5+6*x^4+2*x^3+x^2+2*x+1] 1 30`.
`gen.py` must **skip any line that does not start with a digit** (header,
blank lines), and skip lines whose value field is `?`.

Each subsequent line is one prime. The number of coefficients equals the
**genus g** of the curve (this is the source of truth — from `dump_lpoly` in
`lpdata.c`, which prints `a[0..n-1]` where `n` = number of L-poly coeffs = g):

| genus | line format        | meaning                                            |
|-------|--------------------|----------------------------------------------------|
| 1     | `p,a1`             | a1 = linear coeff of L_p(T) = **-a_p**             |
| 2     | `p,a1,a2`          | a1, a2 = first two coeffs of L_p(T)                |
| 3     | `p,a1,a2,a3`       | (not used by us)                                   |
| any   | `p,?`              | L-poly not computed (bad reduction); **skip it**   |

### Reconstructing the full L-polynomial

smalljac prints only the first `g` coefficients; the rest follow from the
functional equation of the (motivic-weight-1) curve L-poly, which is reciprocal
of degree `2g` with `c_0 = 1`:  `c_{2g-i} = p^{g-i} c_i`.

- **Genus 1** (degree-2 L-poly), `p,a1`:
  `L_p(T) = 1 + a1·T + p·T^2`  → ascending `[1, a1, p]`, with `a1 = -a_p`.
  (This is exactly what `bench/gen.py` already consumes as `(p, t)`, `t = a1`.)

- **Genus 2** (degree-4 L-poly), `p,a1,a2`:
  `L_p(T) = 1 + a1·T + a2·T^2 + a1·p·T^3 + p^2·T^4`
  → ascending **`[1, a1, a2, a1*p, p*p]`**.

This genus-2 reconstruction was checked against the Pari-validated factors in
`bench/gen_g2.py` for curve 169.a and matches exactly (see below).

## Genus-2 curve-spec gotcha (IMPORTANT for objects.yaml)

smalljac parses a genus-2 curve only as a **`y^2 = f(x)` polynomial string in
`x`** (e.g. `"x^6+4*x^5+6*x^4+2*x^3+x^2+2*x+1"`). It does **NOT** accept:
- the LMFDB `[[f],[h]]` pair form (`y^2 + h(x) y = f(x)`), nor
- a bracketed ascending-coefficient list `[a0,a1,...]` for the sextic.
  (The bracket-list form is only for elliptic / number-field inputs, e.g.
  `[1,-1-a,a,0,0] / (a^2-a-1)`.)

So an LMFDB g2c model `[f,h]` must be **completed-the-square** before being
handed to smalljac: `y^2 = 4·f(x) + h(x)^2`. For 169.a, `f=[0,0,0,0,1,1]`,
`h=[1,1,0,1]` gives `4f+h^2 = [1,2,1,2,6,4,1]` = `x^6+4x^5+6x^4+2x^3+x^2+2x+1`,
which is the model used in the worked example below. (Completing the square
does not change the curve's L-function.)

Therefore in `objects.yaml`, the `curve` field for a `genus2` object should be
the smalljac `y^2=f` polynomial string (or the `gen.py` `--backend smalljac`
path must complete the square from the stored `[f,h]` itself).

## Worked example (actual output of our genus-2 build)

Built smalljac v4.1.3 (genus 2) + ff_poly v1.2.7 and ran:

```
$ lpdata g1 "[0,-1,1,-10,-20]" 30 1        # genus 1 (elliptic 11.a1), good-only
# g1_lpdata.txt:
[0,-1,1,-10,-20] 1 30
2,2
3,1
5,-1
7,2
13,-4
17,2
19,0
23,1
29,0

$ lpdata g2s "x^6+4*x^5+6*x^4+2*x^3+x^2+2*x+1" 30 1   # genus 2 (169.a), good-only
# g2s_lpdata.txt:
[x^6+4*x^5+6*x^4+2*x^3+x^2+2*x+1] 1 30
3,2,1
5,0,-7
7,0,7
11,0,11
17,-3,-8
19,6,31
23,-6,13
29,3,-20
```

Reconstructed genus-2 L-polys (ascending `[1,a1,a2,a1*p,p^2]`) vs `gen_g2.py`
Pari `hyperellcharpoly` self-test for 169.a:

| p | smalljac `a1,a2` | reconstructed L-poly | gen_g2 expected | match |
|---|------------------|----------------------|-----------------|-------|
| 3 | 2, 1             | `[1,2,1,6,9]`        | `[1,2,1,6,9]`   | ✓     |
| 5 | 0, -7            | `[1,0,-7,0,25]`      | `[1,0,-7,0,25]` | ✓     |
| 7 | 0, 7             | `[1,0,7,0,49]`       | (not in selftest, consistent) | ✓ |

(Note `lpdata` with good-only skips the bad prime p=13 for 169.a; bad factors
are injected from LMFDB by the driver, never taken from lpdata — per INTERFACES.md §D.)
