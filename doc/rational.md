# The `rational` command-line tool

`examples/rational.c` is a small command-line front end that computes an
L-function from a one-line text spec, without writing any C against
[`glfunc.h`](../include/glfunc.h). It is the quickest way to drive the library
for a motivic L-function of a "rational" object (an elliptic curve, a genus-2
curve, a classical newform, ...) once you have its degree, conductor, motivic
weight, Gamma-shifts, and Euler factors in hand.

For programmatic use against the C API, see the API guide at
[`doc/api.md`](api.md); for symmetric-power L-functions, see
[`doc/sympow.md`](sympow.md).

## Building and invoking

Build the example (it is compiled by `make` / `make executables`):

```sh
./configure
make executables          # or: make build/examples/rational.exe
```

The binary takes exactly **two** arguments, an input file and an output file:

```
./build/examples/rational.exe <input_file> <output_file>
```

It reads `<input_file>` line by line, computing one L-function per non-empty
line, and prints the result of each to **standard output**. The `<output_file>`
argument is currently required by the tool but not yet written to (the
per-line results are reported on stdout); pass any writable path for it.

Each result block printed to stdout contains the analytic rank, the root number
(epsilon, a complex ball), the leading Taylor coefficient at the central point,
and the first zero on the critical line, each as a certified ball (midpoint
`+/-` radius).

## Input-line grammar

Each line has exactly six colon-separated fields:

```
label:degree:conductor:weight:[mus]:[[euler_factors]]
```

A line is tokenised by splitting on `:` into exactly six tokens (a line with a
different number of `:`-separated fields is rejected). The fields are:

| Field | Type | Meaning |
| --- | --- | --- |
| `label` | string | A free-form identifier for the object (echoed back in the output; not otherwise interpreted). It must not contain a `:`. |
| `degree` | integer | The degree `d` of the L-function. This fixes the length of `mus` (`d` entries) and the width of each Euler factor (`d + 1` coefficients). |
| `conductor` | integer | The conductor `N`. Together with the Gamma-shifts it determines `Lfunc_nmax`, the largest prime at which an Euler factor is needed. |
| `weight` | integer | The motivic weight `w`. The tool passes `normalisation = weight/2` to `Lfunc_init` (so an elliptic curve, motivic weight 1, uses `weight = 1`; a classical newform of weight `k` uses motivic weight `k - 1`). |
| `[mus]` | list of `degree` numbers | The Gamma_R shifts `mu_i`, as a bracketed comma-separated list. Each `mu_i + normalisation` must be a non-negative half-integer (otherwise `Lfunc_init` fails). Values may be integers or decimals. |
| `[[euler_factors]]` | list of lists | The local L-polynomial coefficients, one inner list per prime; see below. |

### The `[mus]` field

`[mus]` is a single bracketed list of exactly `degree` numbers, for example
`[0,1]` for an elliptic curve or `[0,0,1,1]` for a genus-2 curve. With
`normalisation = weight/2`, the convention matches `Lfunc_init`: for an elliptic
curve over Q (`weight = 1`, so `normalisation = 0.5`), `mus = [0,1]` gives the
half-integer shifts `[0.5, 1.5]` after normalisation.

### The `[[euler_factors]]` field

`[[euler_factors]]` is a list of inner lists. Each inner list is the coefficient
vector of one local L-polynomial, written in **ascending** order (constant term
first), and must have exactly `degree + 1` entries:

```
[c0, c1, ..., c_degree]   <->   c0 + c1*T + ... + c_degree*T^degree
```

with the convention `c0 = 1`. The inner lists are matched to the primes
`2, 3, 5, 7, 11, ...` **in order**: the first inner list is the factor at
`p = 2`, the second at `p = 3`, and so on. You must supply at least as many
factors as there are primes up to `Lfunc_nmax` (the tool asserts this); for an
elliptic curve, a degree-2 factor is `[1, -a_p, p]` at a good prime `p`.

Coefficients are parsed as 64-bit integers, so this tool is suited to objects
whose local factors have integer coefficients that fit in `int64` (elliptic
curves, genus-2 curves, and newforms with rational coefficients). Symmetric
powers, whose coefficients overflow 64 bits, are assembled differently; see
[`doc/sympow.md`](sympow.md).

The bad primes are included in the list at their prime position, with their
(usually lower-degree) local factor padded to `degree + 1` entries with trailing
zeros. For example, an elliptic curve of conductor 11 has a degree-1 bad factor
at `p = 11`, written as `[1, -a_11, 0]` (a trailing `0` so the inner list still
has `degree + 1 = 3` entries).

## Worked example

The elliptic curve `11.a2` (LMFDB label, equation `y^2 + y = x^3 - x^2`) has
degree 2, conductor 11, motivic weight 1, and `mus = [0,1]`. The file header of
`examples/rational.c` gives its spec with Euler factors at every prime up to 79
(the bad prime 11 appears as `[1,-1,0]`, split multiplicative with `a_11 = 1`):

```
11a2:2:11:1:[0,1]:[[1,2,2],[1,1,3],[1,-1,5],[1,2,7],[1,-1,0],[1,-4,13],[1,2,17],[1,0,19],[1,1,23],[1,0,29],[1,-7,31],[1,-3,37],[1,8,41],[1,6,43],[1,-8,47],[1,6,53],[1,-5,59],[1,-12,61],[1,7,67],[1,3,71],[1,-4,73],[1,10,79]]
```

Save that single line to `input.txt` and run:

```text
$ rm -f g_*                       # scrub any stale G-cache in the CWD first
$ ./build/examples/rational.exe input.txt output.txt
Input: input.txt
Output: output.txt
label = 11a2
degree = 2 conductor = 11 weight = 1 mus = [ 0.00 1.00 ]
Rank = 0
Epsilon = (1.0000000000000000000 + 0j)  +/-  (2.26e-115, 6.72e-58j)
Leading Taylor coeff = 0.25384186085591068434 +/- 1.8368e-57
First zero = 6.3626138947130887014 +/- 8.3524e-53
```

The rank, the root number (`+1`), the leading Taylor coefficient
(`0.25384186...`), and the first zero (`6.36261389...`) agree with the LMFDB
data for this isogeny class (its newform is `11.2.a.a`), and with the `11.a1`
golden values in `test/highdeg/objects.yaml` (`11.a1` and `11.a2` share an
L-function).

## Where the Euler-factor data comes from

The tool consumes Euler factors; it does not compute them (the library has no
point counting or Hecke-eigenvalue machinery of its own). Typical sources:

- **[LMFDB](https://www.lmfdb.org).** Each L-function / object page lists the
  Euler factors (or the `a_p`, from which `[1, -a_p, p]` is formed for an
  elliptic curve). This is the easiest source for a single named object.
- **smalljac / `lpdata`.** Andrew Sutherland's `smalljac` ships an `lpdata`
  utility that emits good-prime local data for elliptic and genus-2 curves in
  bulk; this is what the high-degree regression suite uses (see
  [`test/highdeg/INTERFACES.md`](https://github.com/edgarcosta/lfunctions/blob/main/test/highdeg/INTERFACES.md) and
  [`test/highdeg/SMALLJAC_FORMAT.md`](https://github.com/edgarcosta/lfunctions/blob/main/test/highdeg/SMALLJAC_FORMAT.md)). Bad
  primes are not produced by `lpdata` and must be supplied from a reliable
  source such as the LMFDB.
- **Pari/GP, Sage, Magma.** Any of these can produce `a_p` / local factors for
  the object you have in mind.

Whatever the source, the factors must be in the ascending, `degree + 1`-wide,
prime-ordered form described above, and they must cover every prime up to
`Lfunc_nmax`.

> **Stale `g_*` cache files.** The library caches the gamma-factor product to
> disk as `g_<normalisation>` files in the current directory. Remove any such
> files (`rm -f g_*`) before a run; a stale one is silently reused and can
> corrupt the result.
