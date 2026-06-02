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

2. Change your working directory and run the configure file
```
cd lfunctions
./configure <options>
```

3. If any of these libraries are installed in some other location than the default path `/usr/local`, pass `--with-flint=...`, `--with-primesieve=...`, `--with-gmp=...`, or `--with-mpfr=...` with the correct path to configure (type ./configure --help to show more options).

Note that with exception of `primesieve`, all these libraries are already provided by [SageMath](http://www.sagemath.org/),
and one can do `--with-flint=<SAGE_DIR>/local/`

4. Compile everything by doing
```
make
```

5. You can use it as shared library, the inteface is defined in 
interface is defined in `include/glfunc.h` and the shared library is provided in
`build/liblfun.so` (or `build/liblfun.dylib` on macOS).




