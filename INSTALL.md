# Installing and building lfunctions

This document is the full build and install guide for `lfunctions`, the generic
L-function calculator (shared library `liblfun`). The `README.md` keeps only a
short summary and links here.

The build is driven by [GNU Autoconf](https://www.gnu.org/software/autoconf/):
the build description lives in `configure.ac`, and the generated `./configure`
script is **committed to the repository**, so building from a checkout needs no
autotools. `configure` substitutes its results into `Makefile.in` to produce the
`Makefile`. The committed `configure` was generated with Autoconf 2.72.

## Prerequisites

| Dependency | Minimum version | Notes |
| --- | --- | --- |
| [FLINT](http://flintlib.org/) | 3.0 | Since FLINT 3.0 the Arb ball-arithmetic library is merged into FLINT: the `flint/acb.h`, `flint/arb.h`, ... headers ship in FLINT and the `acb_*` / `arb_*` symbols live in `libflint`, so **no separate Arb install is needed**. (FLINT 2.x with a standalone `libarb` is detected as a fallback, but 3.0+ is the supported configuration.) |
| [primesieve](https://github.com/kimwalisch/primesieve) | any recent | Fast prime generator, used to enumerate the primes at which Euler factors are supplied. |
| [GMP](https://gmplib.org/) | any FLINT-compatible | Pulled in by FLINT. |
| [MPFR](http://www.mpfr.org/) | any FLINT-compatible | Pulled in by FLINT. |
| C / C++ toolchain | C11 (`gnu11`) and C++17 (`c++1z`) | The library is C; the example and test drivers are C++. A recent GCC or Clang works. |

You also need `make`. You do **not** need autotools unless you edit
`configure.ac` (see [Regenerating configure](#regenerating-configure)).

### The SageMath shortcut

[SageMath](http://www.sagemath.org/) already ships FLINT (>= 3.0), GMP, and MPFR.
With a Sage install you only need to add `primesieve` separately, then point
`configure` at Sage's prefix:

```
./configure --with-flint=$SAGE_DIR/local --with-primesieve=/usr/local
```

`$SAGE_DIR` is the root of your Sage installation (the directory containing
`local/`). Because FLINT provides GMP, MPFR, and Arb under that same prefix, the
single `--with-flint` usually pins all three.

## Platform notes

`lfunctions` builds on Linux and macOS. `configure` adapts to the host OS:

- **Shared-library extension.** The library is `build/liblfun.so` on Linux and
  `build/liblfun.dylib` on macOS. The Makefile variable `SHLIB_EXT` (`so` or
  `dylib`) is set by `configure` from `uname -s`; the rest of this document
  writes `.so` for brevity.
- **The `.exe` suffix is used on every platform, including Linux and macOS.**
  Compiled tests live in `build/test/*.exe` and compiled examples in
  `build/examples/*.exe`. This is a naming convention of this project's
  Makefile, not a Windows artifact.
- **rpath.** When a dependency is supplied via `--with-<name>=DIR`, that
  library directory is baked into the binaries as an rpath, so the examples and
  tests run without setting `LD_LIBRARY_PATH` (Linux) or
  `DYLD_LIBRARY_PATH` (macOS). On Linux the library additionally links `-lrt`;
  on macOS the install name is `@rpath/liblfun.dylib`.

## 1. Get the source

```
git clone https://github.com/edgarcosta/lfunctions.git
cd lfunctions
```

## 2. Configure

```
./configure [options]
```

`configure` locates the dependencies, then writes `Makefile` from `Makefile.in`
and prints a configuration report (host OS, chosen compilers, the resolved
header and library paths for each dependency, and any search directories added).
Re-run `./configure` after editing `Makefile.in` or changing any flag.

### Finding the dependencies

`configure` looks for each dependency in this order:

1. Any prefix you pass explicitly with `--with-<name>=DIR` (searched first).
2. `pkg-config`, if it is installed: `configure` runs `pkg-config` for `flint`
   and `primesieve` and adds their reported include and library directories
   automatically.
3. The compiler and linker default search paths (typically `/usr` and
   `/usr/local`).

If a library lives outside the default search path and `pkg-config` does not
know about it, pass the matching `--with-<name>=DIR`, where `DIR` is the install
prefix (the directory containing `lib/` and `include/`).

### Configure options

These are the package-specific options (run `./configure --help` for the full
Autoconf-generated list, including the standard `--prefix`, `--bindir`, and
environment variables such as `CC`, `CFLAGS`, `CPPFLAGS`, `LDFLAGS`, `LIBS`,
`CXX`, `CXXFLAGS`):

| Option | Effect |
| --- | --- |
| `--with-flint=DIR` | Use the FLINT (and bundled Arb) installed under prefix `DIR`. |
| `--with-primesieve=DIR` | Use the primesieve installed under prefix `DIR`. |
| `--with-gmp=DIR` | Use the GMP installed under prefix `DIR`. |
| `--with-mpfr=DIR` | Use the MPFR installed under prefix `DIR`. |
| `--disable-assert` | Compile with `-DNDEBUG`, which turns `assert`s **off**. Asserts are **on** by default, and the test suite uses them as its oracle, so leave them on when running the tests. |
| `--disable-gdb` | Omit the `-g` debug symbols that are added by default. |
| `--enable-sanitize` | Build with `-fsanitize=address,undefined -fno-omit-frame-pointer` (ASan / UBSan). Intended for CI and debugging; off by default. The sanitizer flags are appended to both the compile and link steps, so `libasan` must be present at link time. |

The default compiler flags are
`-pedantic -Wall -Wextra -O2 -funroll-loops -fPIC` (plus `-g` unless
`--disable-gdb`, plus `-DNDEBUG` if `--disable-assert`), with `-std=gnu11` for C
and `-std=c++1z` for C++. Setting `CFLAGS` or `CXXFLAGS` on the `configure` line
**replaces** these defaults (the `-fPIC` needed for the shared link is part of
the default set, so if you override `CFLAGS` keep `-fPIC`).

## 3. Build

All build artifacts live under `build/`. The make targets are:

- `make` (or `make all`): build the object files, the shared library, the
  tests, and the examples.
- `make lib`: build only the shared library `build/liblfun.so` (`.dylib` on
  macOS).
- `make test`: *compile* the tests into `build/test/*.exe` (it does not run
  them).
- `make executables`: compile the examples into `build/examples/*.exe`.
- `make clean`: remove the `build/` directory (`rm -rf build`).

`make` runs in parallel by default: `configure` sets the `JOBS` variable to the
CPU count plus one and the Makefile adds `-j$(JOBS)`.

To run a single test or example, execute its binary directly, for example
`./build/test/dir_test.exe` or `./build/examples/ec_37.a1.exe`.

## 4. Run the tests

### `make check` (the curated regression suite)

```
make check
```

`make check` builds everything, then runs every compiled `test/` and `examples/`
binary (minus a small deny-list of binaries that require command-line
arguments). The **exit code is the oracle**: it exits non-zero if any binary
fails, and the binaries assert against certified values. Because the asserts are
the oracle, run the suite with assertions **on** (the default); a build
configured with `--disable-assert` will compile the asserts out and the tests
will silently pass.

`make check` scrubs stale `g_<...>` cache files before the suite and after each
binary (see [Troubleshooting](#troubleshooting)).

### `make check-highdeg` (the high-degree suite)

```
make check-highdeg
```

This is a separate suite that certifies degree 2 through 9 L-functions
(elliptic curves, their symmetric powers, genus-2 curves, classical newforms)
against LMFDB and Pari golden values. By default it runs purely from committed
base fixtures, so it needs no external toolchain. Regenerating those fixtures, or
running without them, needs `smalljac`'s `lpdata` (and, for some objects, Sage);
in particular the genus-2 objects need a genus-2 `lpdata` build passed as
`LPDATA=.../lpdata2`. Its knobs (`FIXTURES`, `LABEL`, `BACKEND`, `LPDATA`) and
the fixture-regeneration target (`make highdeg-data`) are documented in
[`test/highdeg/INTERFACES.md`](test/highdeg/INTERFACES.md) and the README's
Build & Test section.

## 5. Use the library

The public interface is the single header [`include/glfunc.h`](include/glfunc.h);
the shared library is `build/liblfun.so` (`.dylib` on macOS). To compile and
link your own program against it:

```
cc myprog.c \
  -I/path/to/lfunctions/include \
  $(pkg-config --cflags flint primesieve) \
  -L/path/to/lfunctions/build -llfun \
  $(pkg-config --libs flint primesieve) \
  -o myprog
```

At run time the loader must find `liblfun.so` and the dependency libraries; the
simplest options are to set `LD_LIBRARY_PATH` (Linux) or `DYLD_LIBRARY_PATH`
(macOS) to include `build/`, or to pass an rpath at link time
(`-Wl,-rpath,/path/to/lfunctions/build`).

`examples/rational.c` is a ready-made command-line front end that reads an
L-function spec on one line and computes it without writing any C; see
[`doc/rational.md`](doc/rational.md). The C examples (for instance
`examples/rational.c` and `test/dir_test.c`) use only `glfunc.h`; the C++
examples additionally use FLINT's C++ bindings.

## Troubleshooting

- **Header / library prefix-mismatch warning.** If a dependency's header and
  library are found under different install prefixes, `configure` prints a
  non-fatal warning of the form

  ```
  configure: WARNING: flint header and library resolve to different prefixes --
  configure: WARNING:   you may compile against one flint and link another:
  ```

  This means the build would compile against one copy of the dependency and link
  another, a silent ABI mismatch. Pass `--with-<name>=DIR` to pin both the
  header and the library to the same prefix. (This commonly happens when a
  system copy and a Sage or local copy are both installed.)

- **Stale `g_<...>` cache files poison results.** The G-function (the gamma
  factor product) is cached to disk in the current working directory as files
  named `g_<normalisation>` (for example `g_0.5`, `g_1.5`). A stale or
  mismatched cache file sitting in the CWD is silently reused and can corrupt a
  computation. `make check` and `make check-highdeg` scrub `g_[0-9]*` for you,
  but for ad-hoc runs remove them first:

  ```
  rm -f g_*
  ```

  For hermetic runs, point `cache_dir` (in `Lparams_t`, via
  `Lfunc_init_advanced`) at a clean, per-run directory.

- **`make` does not pick up a `Makefile.in` or flag change.** The `Makefile` is
  generated. Re-run `./configure` (with the same options) after editing
  `Makefile.in` or changing any configure flag, then `make` again.

- **A dependency is installed but not found.** Pass its prefix with
  `--with-<name>=DIR`. The configuration report at the end of `./configure`
  prints the exact header and library paths the toolchain selected for each
  dependency, which is the quickest way to see what was picked up.

## Regenerating configure

You only need this if you edit `configure.ac`; the generated `./configure` is
committed, so ordinary builds never regenerate it. Maintainers can rebuild it
with:

```
make configure
```

which runs `autoreconf -i` and therefore requires `autoconf` to be installed.
