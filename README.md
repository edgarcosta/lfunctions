# Generic L-function calculator

A generic L-function calculator for motivic L-functions, originally written by Dave J. Platt.

For details of the methods we refer to:
 - *Computing motivic L-functions*, (under preparation) by Jonathan W. Bober, Andrew R. Booker, Edgar Costa, Min Lee, David J. Platt, and Andrew V. Sutherland.


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
 - `make clean`: remove the `build/` directory (`rm -rf build`).

To run a single test or example, execute its binary directly, e.g. `./build/test/dir_test.exe` or `./build/examples/ec_37.a1.exe`.

See [INSTALL.md](INSTALL.md) for the configure options, the SageMath shortcut, and troubleshooting (including scrubbing stale `g_*` cache files).

