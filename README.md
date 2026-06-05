# Generic L-function calculator

A generic L-function calculator for motivic L-functions, originally written by Dave J. Platt.

For details of the methods we refer to:
 - *Computing motivic L-functions*, (under preparation) by Jonathan W. Bober, Andrew R. Booker, Edgar Costa, Min Lee, David J. Platt, and Andrew V. Sutherland.


## Documentation

The API reference and guides are built with [Sphinx](https://www.sphinx-doc.org/)
from the sources under `doc/`. The per-function reference is generated from the
comments in `include/glfunc.h` (via Doxygen and Breathe), so the header stays the
single source of truth. The docs build is self-contained under `doc/` and is
independent of the library build.

Build the HTML locally, from the repository root:

```sh
make -C doc html
```

This writes the site to `doc/_build/html/index.html`. It needs `sphinx-build`
(with the `myst-parser` and `breathe` packages) and the `doxygen` binary. A
convenient setup is a virtual environment under `doc/`:

```sh
python3 -m venv doc/.venv
doc/.venv/bin/pip install -r doc/requirements.txt
sudo apt-get install doxygen        # or: brew install doxygen
make -C doc html
```

`make -C doc html` (equivalently `cd doc && make html`) runs Doxygen over
`include/glfunc.h` and then Sphinx. `doc/Makefile` is a standard Sphinx
makefile, so its other targets work too (`make -C doc clean`, `make -C doc help`).


## Dependencies

The primary dependencies are:

 - [FLINT: Fast Library for Number Theory](http://flintlib.org/), version ≥ 3.0. Since FLINT 3.0 the Arb ball-arithmetic library is merged into FLINT (the `flint/acb.h`, `flint/arb.h`, ... headers live in FLINT, and the `acb_*`/`arb_*` symbols are in `libflint`), so **no separate Arb install is needed**.
 - [primesieve - Fast prime number generator](https://github.com/kimwalisch/primesieve)
 - [GMP](https://gmplib.org/) and [MPFR](http://www.mpfr.org/), which FLINT pulls in.

[SageMath](http://www.sagemath.org/) already ships FLINT (≥ 3.0), GMP, and MPFR; with Sage you only need to install `primesieve` separately.

## Install

The build is driven by [GNU Autoconf](https://www.gnu.org/software/autoconf/); the generated `./configure` is committed, so a checkout needs no autotools. The short version is:

```
git clone https://github.com/edgarcosta/lfunctions.git
cd lfunctions
./configure                 # add --with-flint=DIR etc. if deps are not on the default path
make
```

The public interface is `include/glfunc.h`; the shared library is built at `build/liblfun.so` (`build/liblfun.dylib` on macOS).

For the full guide, see **[INSTALL.md](INSTALL.md)**: prerequisites and versions, platform notes, every `./configure` option, the SageMath shortcut, linking against `liblfun`, and troubleshooting.

## Build & Test

`configure` substitutes its results into `Makefile.in` to produce the `Makefile` (the `Makefile` is generated; do not edit it). Re-run `./configure` after editing `Makefile.in` or changing any configure flag. All build artifacts live under `build/`, and the binaries use a `.exe` suffix on every platform, including Linux and macOS.

The main targets are:

 - `make` (or `make all`): build the object files, the shared library, the tests, and the examples.
 - `make lib`: build only the shared library `build/liblfun.so` (`.dylib` on macOS).
 - `make test` / `make executables`: *compile* the tests / examples into `build/test/*.exe` and `build/examples/*.exe` (they do not run them).
 - `make check`: build everything, then run the curated self-checking regression suite; it exits non-zero if any binary fails (the exit code is the oracle). It also scrubs stale `g_<...>` cache files before the suite and after each binary.
 - `make check-highdeg`: run the high-degree regression suite (certified degree 2 to 9 L-functions: elliptic curves, their symmetric powers, genus-2 curves, classical newforms) against LMFDB and Pari golden values. By default it runs purely from committed base fixtures, so it needs no external toolchain. See below.
 - `make clean`: remove the `build/` directory (`rm -rf build`).

To run a single test or example, execute its binary directly, e.g. `./build/test/dir_test.exe` or `./build/examples/ec_37.a1.exe`.

### High-degree regression suite

`make check-highdeg` certifies the degree 2 to 9 objects listed in `test/highdeg/objects.yaml`. By default it runs from the committed base fixtures in `test/highdeg/fixtures/`, so it needs no `smalljac` or Sage; regenerating those fixtures, or running without them, needs `smalljac`'s `lpdata` (and, for the genus-2 objects, a genus-2 `lpdata` build passed as `LPDATA=.../lpdata2`). Its knobs (`FIXTURES`, `LABEL`, `BACKEND`, `LPDATA`), the `make highdeg-data` fixture-regeneration target, and how to add an object to `objects.yaml` are documented in **[test/highdeg/INTERFACES.md](test/highdeg/INTERFACES.md)**.

See [INSTALL.md](INSTALL.md) for the configure options, the SageMath shortcut, and troubleshooting (including scrubbing stale `g_*` cache files).

## Usage

The quickest way to compute an L-function without writing C is the `rational` command-line tool (`examples/rational.c`), which reads a one-line spec `label:degree:conductor:weight:[mus]:[[euler_factors]]` and prints the rank, root number, leading Taylor coefficient, and first zero. See [doc/rational.md](doc/rational.md) for the input grammar and a worked example.

For C callers, the public API accepts Euler factors by callback, one-prime
push, or consecutive-prime array, and also accepts raw Dirichlet coefficient
arrays with an explicit `ALGEBRAIC_NORM` / `ANALYTIC_NORM` selector. See
[doc/api.md](doc/api.md) for the supply contracts and edge-case warnings.

### Symmetric powers

The library can compute symmetric-power L-functions `Sym^k(E)` of an elliptic curve. The good-prime Euler factors are formed from the curve's `a_p` by the helper `sym_power_lpoly` (`include/sym_power.h`); a worked example (Sym^2 and Sym^3 of `11.a1`, with certified assertions) is [`examples/ec_sym.cpp`](examples/ec_sym.cpp). See [doc/sympow.md](doc/sympow.md) for the assembly conventions (degree, normalisation, `mus`, `self_dual`) and the bad-factor caveats.
