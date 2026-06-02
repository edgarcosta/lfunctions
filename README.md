# Generic L-function calculator

A generic L-function calculator for motivic L-functions, originally written by Dave J. Platt.

For details of the methods we refer to:
 - *Computing motivic L-functions*, (under preparation) by Jonathan W. Bober, Andrew R. Booker, Edgar Costa, Min Lee, David J. Platt, and Andrew V. Sutherland.


## Dependencies

The primary dependencies are:

 - [FLINT: Fast Library for Number Theory](http://flintlib.org/) — version ≥ 3.0. Since FLINT 3.0 the Arb ball-arithmetic library is merged into FLINT (the `flint/acb.h`, `flint/arb.h`, … headers live in FLINT, and the `acb_*`/`arb_*` symbols are in `libflint`), so **no separate Arb install is needed**.
 - [primesieve - Fast prime number generator](https://github.com/kimwalisch/primesieve)

FLINT in turn pulls in:

 - [GMP: GNU Multiple Precision Arithmetic Library](https://gmplib.org/)
 - [MPFR: GNU Multiple Precision Floating-Point Reliably](http://www.mpfr.org/)

[SageMath](http://www.sagemath.org/) already ships FLINT (≥ 3.0), GMP, and MPFR; with Sage you only need to install `primesieve` separately. 

## Installation

1. Download `lfunctions`
```
git clone https://github.com/edgarcosta/lfunctions.git
```

2. Change your working directory and run the configure script
```
cd lfunctions
./configure <options>
```

3. If any of these libraries are installed in some other location than the default path `/usr/local`, pass `--with-flint=DIR`, `--with-primesieve=DIR`, `--with-gmp=DIR`, or `--with-mpfr=DIR` with the correct prefix (run `./configure --help` to see all options). When `pkg-config` is available, `configure` will also pick up FLINT and primesieve automatically.

Note that, with the exception of `primesieve`, all these libraries are already provided by [SageMath](http://www.sagemath.org/),
so with a Sage install you can do `--with-flint=<SAGE_DIR>/local/`.

4. Compile everything by doing
```
make
```

5. You can use it as a shared library: the interface is defined in `include/glfunc.h`, and the shared library is built at
`build/liblfun.so` (or `build/liblfun.dylib` on macOS).

## Build & Test

The build is driven by [Autoconf](https://www.gnu.org/software/autoconf/). The build description lives in `configure.ac`, and the generated `./configure` script is **committed to the repository**, so building from a checkout needs no autotools — just run `./configure` as above. `configure` substitutes its results into `Makefile.in` to produce the `Makefile` (the `Makefile` itself is generated; do not edit it). Re-run `./configure` after editing `Makefile.in` or changing any configure flag.

Useful `./configure` flags:

 - `--with-flint=DIR`, `--with-primesieve=DIR`, `--with-gmp=DIR`, `--with-mpfr=DIR` — prefixes for dependencies installed outside the default search path.
 - `--disable-assert` — compile with `-DNDEBUG` (turns `assert`s off; they are **on** by default).
 - `--disable-gdb` — omit `-g` debug symbols (added by default).
 - `--enable-sanitize` — build with `-fsanitize=address,undefined` (ASan/UBSan; intended for CI/debugging, off by default).

Maintainers who edit `configure.ac` can regenerate the committed `./configure` with `make configure` (this runs `autoreconf -i` and so requires `autoconf` to be installed).

All build artifacts live under `build/` (the binaries use a `.exe` suffix on every platform, including Linux). The make targets are:

 - `make` (or `make all`) — build the object files, the shared library, the tests, and the examples.
 - `make lib` — build only the shared library `build/liblfun.so` (`.dylib` on macOS).
 - `make test` — *compile* the tests into `build/test/*.exe` (it does not run them).
 - `make executables` — compile the examples into `build/examples/*.exe`.
 - `make check` — build everything, then run the curated self-checking regression suite; it exits non-zero if any binary fails (the exit code is the oracle). It also scrubs stale `g_<…>` cache files before the suite and after each binary.
 - `make clean` — remove the `build/` directory (`rm -rf build`).

To run a single test or example, execute its binary directly, e.g. `./build/test/dir_test.exe` or `./build/examples/ec_37.a1.exe`.

