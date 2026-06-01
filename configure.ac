#                                               -*- Autoconf -*-
# Process with `autoreconf -i` (or `make configure`) to produce ./configure.
#  Copyright (C) 2016-2026 Edgar Costa
#  See the file LICENSE for license details.

AC_PREREQ([2.69])
AC_INIT([lfunctions], [0.1])
AC_CONFIG_SRCDIR([include/glfunc.h])

# Record whether the user supplied CFLAGS/CXXFLAGS *before* AC_PROG_CC touches them.
# Use ${VAR:+yes} (not ${VAR+yes}) so a set-but-empty value is treated as "unset" and still
# gets the defaults -- matching the old configure's [ -z "$CFLAGS" ] test, so `CFLAGS= ./configure`
# does not silently drop -fPIC (which would break the -shared link).
lfun_user_cflags_set=${CFLAGS:+yes}
lfun_user_cxxflags_set=${CXXFLAGS:+yes}
lfun_user_cflags=$CFLAGS
lfun_user_cxxflags=$CXXFLAGS

AC_PROG_CC
AC_PROG_CXX

# ---- options -------------------------------------------------------------
AC_ARG_ENABLE([assert],
  [AS_HELP_STRING([--disable-assert], [compile with -DNDEBUG (turns asserts off)])],
  [], [enable_assert=yes])
AC_ARG_ENABLE([gdb],
  [AS_HELP_STRING([--disable-gdb], [do not add -g debug symbols])],
  [], [enable_gdb=yes])

AC_ARG_WITH([flint],      [AS_HELP_STRING([--with-flint=DIR],      [location of FLINT (>=3.0; provides Arb)])])
AC_ARG_WITH([gmp],        [AS_HELP_STRING([--with-gmp=DIR],        [location of GMP])])
AC_ARG_WITH([mpfr],       [AS_HELP_STRING([--with-mpfr=DIR],       [location of MPFR])])
AC_ARG_WITH([primesieve], [AS_HELP_STRING([--with-primesieve=DIR], [location of primesieve])])

# ---- collect include/lib dirs from --with and pkg-config -----------------
INC_DIRS=
LIB_DIRS=

# add_dep_dir NAME DIR : find lib subdir (installed or in-tree) + include dir
add_dep_dir () {
  lfun_name=$1
  lfun_dir=$2
  if   test -d "$lfun_dir/lib";       then lfun_libdir=$lfun_dir/lib;       lfun_incdir=$lfun_dir/include
  elif test -d "$lfun_dir/lib64";     then lfun_libdir=$lfun_dir/lib64;     lfun_incdir=$lfun_dir/include
  elif test -d "$lfun_dir/.libs";     then lfun_libdir=$lfun_dir/.libs;     lfun_incdir=$lfun_dir
  elif test -d "$lfun_dir/src/.libs"; then lfun_libdir=$lfun_dir/src/.libs; lfun_incdir=$lfun_dir/src
  else echo "configure: error: invalid $lfun_name directory: $lfun_dir" >&2; exit 1
  fi
  LIB_DIRS="$LIB_DIRS $lfun_libdir"
  INC_DIRS="$INC_DIRS $lfun_incdir"
}

for pair in "FLINT:$with_flint" "GMP:$with_gmp" "MPFR:$with_mpfr" "PRIMESIEVE:$with_primesieve"; do
  lfun_n=${pair%%:*}
  lfun_v=${pair#*:}
  case "$lfun_v" in
    ""|yes|no) ;;                       # absent / bare / --without : rely on default search
    *) add_dep_dir "$lfun_n" "$lfun_v" ;;
  esac
done

AC_CHECK_PROG([PKG_CONFIG], [pkg-config], [pkg-config])
if test -n "$PKG_CONFIG"; then
  for m in flint primesieve; do
    if $PKG_CONFIG --exists "$m" 2>/dev/null; then
      for i in `$PKG_CONFIG --cflags-only-I "$m"`; do INC_DIRS="$INC_DIRS ${i#-I}"; done
      for l in `$PKG_CONFIG --libs-only-L "$m"`;  do LIB_DIRS="$LIB_DIRS ${l#-L}"; done
    fi
  done
fi

# ---- dependency detection (seed search with the collected dirs) ----------
lfun_inc_flags=; for d in $INC_DIRS; do lfun_inc_flags="$lfun_inc_flags -I$d"; done
lfun_lib_flags=; for d in $LIB_DIRS; do lfun_lib_flags="$lfun_lib_flags -L$d"; done

save_CPPFLAGS=$CPPFLAGS; save_LDFLAGS=$LDFLAGS; save_LIBS=$LIBS
CPPFLAGS="$CPPFLAGS $lfun_inc_flags"
LDFLAGS="$LDFLAGS $lfun_lib_flags"

AC_LANG_PUSH([C])
AC_CHECK_HEADER([gmp.h], [], [AC_MSG_ERROR([gmp.h not found; install GMP or pass --with-gmp=DIR])])
AC_CHECK_LIB([gmp], [__gmpz_init], [], [AC_MSG_ERROR([libgmp not found; install GMP or pass --with-gmp=DIR])])
AC_CHECK_HEADER([mpfr.h], [], [AC_MSG_ERROR([mpfr.h not found; install MPFR or pass --with-mpfr=DIR])])
AC_CHECK_LIB([mpfr], [mpfr_init], [], [AC_MSG_ERROR([libmpfr not found; install MPFR or pass --with-mpfr=DIR])], [-lgmp])
AC_CHECK_HEADER([primesieve.h], [], [AC_MSG_ERROR([primesieve.h not found; install primesieve or pass --with-primesieve=DIR])])
AC_CHECK_LIB([primesieve], [primesieve_generate_primes], [], [AC_MSG_ERROR([libprimesieve not found; install primesieve or pass --with-primesieve=DIR])])
AC_CHECK_HEADER([flint/flint.h], [], [AC_MSG_ERROR([flint/flint.h not found; install FLINT >=3.0 or pass --with-flint=DIR])])
AC_CHECK_HEADER([flint/acb.h], [], [AC_MSG_ERROR([flint/acb.h not found; need FLINT >=3.0 (Arb merged)])])
AC_CHECK_LIB([flint], [fmpz_init], [], [AC_MSG_ERROR([libflint not found; install FLINT >=3.0 or pass --with-flint=DIR])], [-lmpfr -lgmp])
# acb_init lives in libflint (FLINT>=3) or libarb (FLINT 2.x). Never blindly -larb.
lfun_arb_lib=
AC_SEARCH_LIBS([acb_init], [flint arb],
  [AS_IF([test "x$ac_cv_search_acb_init" = "x-larb"], [lfun_arb_lib=arb])],
  [AC_MSG_ERROR([acb_init not found in libflint or libarb; need FLINT >=3.0 or Arb])],
  [-lflint -lmpfr -lgmp])
AC_LANG_POP([C])

CPPFLAGS=$save_CPPFLAGS; LDFLAGS=$save_LDFLAGS; LIBS=$save_LIBS

# ---- compose INCS --------------------------------------------------------
# own headers -I; dependency headers -isystem, deduped, skipping standard dirs
# (forcing -isystem/usr/include ahead of toolchain dirs breaks libstdc++'s
#  #include_next <stdlib.h>).
INCS="-I'\$(CURDIR)/include'"
lfun_seen=
for d in $INC_DIRS; do
  case " $lfun_seen " in *" $d "*) continue ;; esac
  lfun_seen="$lfun_seen $d"
  case "$d" in /usr/include|/usr/local/include) continue ;; esac
  INCS="$INCS -isystem$d"
done

# ---- compose LIBS --------------------------------------------------------
LLIBS="-L'\$(CURDIR)'"
for d in $LIB_DIRS; do LLIBS="$LLIBS -L$d"; done
lLIBS="-lmpfr -lprimesieve -lgmp"
AS_IF([test -n "$lfun_arb_lib"], [lLIBS="$lLIBS -l$lfun_arb_lib"])
lLIBS="$lLIBS -lflint -lm"

# ---- shared library + rpath, per host OS ---------------------------------
RPATH_FLAGS=
for d in $LIB_DIRS; do RPATH_FLAGS="$RPATH_FLAGS -Wl,-rpath,$d"; done
lfun_os=`uname -s`
case $lfun_os in
  Darwin*) SHLIB_EXT=dylib; SHLIB_LDFLAGS="-dynamiclib -install_name @rpath/liblfun.dylib"; PLATFORM_LIBS= ;;
  *)       SHLIB_EXT=so;    SHLIB_LDFLAGS="-shared -Wl,-soname,liblfun.so";                 PLATFORM_LIBS="-lrt" ;;
esac
LIBS="$LLIBS $lLIBS $PLATFORM_LIBS"

# ---- compose CFLAGS / CXXFLAGS (user value replaces defaults; faithful) --
lfun_base="-pedantic -Wall -Wextra -O2 -funroll-loops -fPIC"
AS_IF([test "x$enable_gdb"    != xno], [lfun_base="$lfun_base -g"])
AS_IF([test "x$enable_assert" =  xno], [lfun_base="$lfun_base -DNDEBUG"])
AS_IF([test "x$lfun_user_cflags_set"   = xyes], [LFUN_CFLAGS=$lfun_user_cflags],     [LFUN_CFLAGS="$lfun_base -std=gnu11"])
AS_IF([test "x$lfun_user_cxxflags_set" = xyes], [LFUN_CXXFLAGS=$lfun_user_cxxflags], [LFUN_CXXFLAGS="$lfun_base -std=c++1z"])

# ---- parallel jobs = ncpu + 1 (replaces the python3 call) ----------------
lfun_ncpu=`getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1`
JOBS=`expr $lfun_ncpu + 1`

AC_SUBST([INCS])
AC_SUBST([LIBS])
AC_SUBST([LFUN_CFLAGS])
AC_SUBST([LFUN_CXXFLAGS])
AC_SUBST([JOBS])
AC_SUBST([SHLIB_EXT])
AC_SUBST([SHLIB_LDFLAGS])
AC_SUBST([RPATH_FLAGS])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

AC_MSG_NOTICE([lfunctions configured:
    os         : $lfun_os
    CC / CXX   : $CC / $CXX
    CFLAGS     : $LFUN_CFLAGS
    CXXFLAGS   : $LFUN_CXXFLAGS
    shared lib : liblfun.$SHLIB_EXT  ($SHLIB_LDFLAGS)])
