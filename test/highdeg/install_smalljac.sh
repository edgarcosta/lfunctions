#!/usr/bin/env bash
#
# install_smalljac.sh — build Andrew Sutherland's `lpdata` (from smalljac),
# configured for genus 2, for the high-degree L-function regression suite.
#
# smalljac v4.1.3 with SMALLJAC_GENUS=2 handles genus 1 AND genus 2 curves
# simultaneously (its smalljac.h: "genus 1 and 2 ... work when SMALLJAC_GENUS
# is set to 2"), which is exactly what this suite needs (genus 1 for elliptic
# curves + symmetric powers, genus 2 for the degree-4 curves). We do NOT need
# genus 3, so we avoid smalljac v5.0 — v5.0 pulls in the `hwlpoly` library
# (-lhwlpoly), which is needed only for genus 3 and is not publicly downloadable.
# v4.1.3 depends only on ff_poly + GMP.
#
# The smalljac / ff_poly sources are NOT on GitHub; they are tarballs on
# https://math.mit.edu/~drew/ . We pin exact tarball filenames (= exact versions)
# below; there are no git tags/commits to pin to.
#
# Output: builds (and optionally installs) the `lpdata` binary and prints its
# absolute path as the FINAL line of stdout, e.g.
#     LPDATA=/path/to/smalljac/smalljac_v4.1.3/lpdata
# A consumer can capture it with:  LPDATA=$(test/highdeg/install_smalljac.sh | tail -1 | cut -d= -f2)
#
# Idempotent: re-running skips downloads/builds that are already present unless
# FORCE=1 is set. CI-friendly: set -euo pipefail, no interactive prompts.
#
# Environment knobs:
#   SMALLJAC_DIR   working dir for sources/builds   (default: <this-dir>/smalljac)
#   PREFIX         where ff_poly headers+lib install (default: $SMALLJAC_DIR/prefix)
#   INSTALL_BIN    if set (e.g. /usr/local/bin), also copy lpdata there (needs write perms)
#   JOBS           parallel make jobs               (default: nproc)
#   FORCE          if 1, wipe and rebuild from scratch
#   FF_POLY_TAR / SMALLJAC_TAR  override the tarball filenames (advanced)

set -euo pipefail

# --- configuration: pinned versions ---------------------------------------
DREW_URL="https://math.mit.edu/~drew"
FF_POLY_TAR="${FF_POLY_TAR:-ff_poly_v1.2.7.tar}"     # ff_poly 1.2.7 (portable CFLAGS, no zn_poly)
SMALLJAC_TAR="${SMALLJAC_TAR:-smalljac_v4.1.3.tar}"  # smalljac 4.1.3 (genus 2, deps: ff_poly + GMP)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SMALLJAC_DIR="${SMALLJAC_DIR:-$SCRIPT_DIR/smalljac}"
PREFIX="${PREFIX:-$SMALLJAC_DIR/prefix}"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
FORCE="${FORCE:-0}"

# Portable compile flags. The upstream makefiles hardcode -O3 with aggressive,
# sometimes non-portable flags (ff_poly v2.0.0 even bakes in -march=icelake-server);
# v1.2.7 is portable, but we still pin our own flags so the build is reproducible
# on any x86-64 CI runner.
PORTABLE_CFLAGS="-O2 -fomit-frame-pointer -funroll-loops -m64 -std=gnu99"

# smalljac v4.1.3 compiles a few functions before their declarations; with modern
# gcc that is an error by default, so allow it as a warning.
SMALLJAC_CFLAGS="$PORTABLE_CFLAGS -Wno-implicit-function-declaration"

log()  { printf '>>> %s\n' "$*" >&2; }
have() { command -v "$1" >/dev/null 2>&1; }

# --- sanity: required tools ------------------------------------------------
for tool in gcc make ar ranlib tar; do
  have "$tool" || { echo "ERROR: required tool '$tool' not found in PATH" >&2; exit 1; }
done
# need a downloader
if have curl;  then DL() { curl -fsSL -o "$2" "$1"; }
elif have wget; then DL() { wget -q -O "$2" "$1"; }
else echo "ERROR: need curl or wget to download tarballs" >&2; exit 1
fi

if [ "$FORCE" = "1" ]; then
  log "FORCE=1: removing $SMALLJAC_DIR"
  rm -rf "$SMALLJAC_DIR"
fi

mkdir -p "$SMALLJAC_DIR" "$PREFIX/include" "$PREFIX/lib"
cd "$SMALLJAC_DIR"

FF_POLY_SRC="$SMALLJAC_DIR/${FF_POLY_TAR%.tar}"      # e.g. .../ff_poly_v1.2.7
SMALLJAC_SRC="$SMALLJAC_DIR/${SMALLJAC_TAR%.tar}"    # e.g. .../smalljac_v4.1.3
LPDATA_BIN="$SMALLJAC_SRC/lpdata"

# --- short-circuit if already built ----------------------------------------
if [ "$FORCE" != "1" ] && [ -x "$LPDATA_BIN" ]; then
  log "lpdata already built at $LPDATA_BIN (set FORCE=1 to rebuild)"
else
  # --- download --------------------------------------------------------------
  for tarname in "$FF_POLY_TAR" "$SMALLJAC_TAR"; do
    if [ ! -f "$SMALLJAC_DIR/$tarname" ]; then
      log "downloading $DREW_URL/$tarname"
      DL "$DREW_URL/$tarname" "$SMALLJAC_DIR/$tarname"
    else
      log "already have $tarname"
    fi
  done

  # --- extract ---------------------------------------------------------------
  [ -d "$FF_POLY_SRC" ]  || { log "extracting $FF_POLY_TAR";  tar xf "$SMALLJAC_DIR/$FF_POLY_TAR"  -C "$SMALLJAC_DIR"; }
  [ -d "$SMALLJAC_SRC" ] || { log "extracting $SMALLJAC_TAR"; tar xf "$SMALLJAC_DIR/$SMALLJAC_TAR" -C "$SMALLJAC_DIR"; }

  # NOTE on stdout discipline: build commands below send their normal output to
  # stderr (>&2) so this script's *stdout* carries ONLY the final `LPDATA=<path>`
  # line. That keeps `LPDATA=$(install_smalljac.sh)` robust regardless of make
  # chatter. Progress still streams live to the terminal / CI log via stderr.

  # --- build + install ff_poly into our private PREFIX -----------------------
  log "building ff_poly ($FF_POLY_TAR)"
  make -C "$FF_POLY_SRC" clean >/dev/null 2>&1 || true
  make -C "$FF_POLY_SRC" -j"$JOBS" CFLAGS="$PORTABLE_CFLAGS" >&2
  # `make install` copies ff_poly.h, the ff_poly/ header dir, and libff_poly.a
  # under $INSTALL_ROOT/{include,lib}.
  make -C "$FF_POLY_SRC" install INSTALL_ROOT="$PREFIX" >&2
  [ -f "$PREFIX/lib/libff_poly.a" ] || { echo "ERROR: ff_poly did not install libff_poly.a" >&2; exit 1; }

  # --- verify smalljac is configured for genus 2 -----------------------------
  # v4.1.3 ships with SMALLJAC_GENUS=2 already; assert it (and fix it if a future
  # tarball ever differs) so the build is unambiguously genus-2.
  if grep -Eq '^[[:space:]]*#define[[:space:]]+SMALLJAC_GENUS[[:space:]]+2\b' "$SMALLJAC_SRC/smalljac.h"; then
    log "smalljac.h has SMALLJAC_GENUS = 2 (genus 1 + 2 supported)"
  else
    log "WARNING: SMALLJAC_GENUS != 2 in smalljac.h; forcing it to 2"
    # portable in-place edit without GNU-sed assumptions
    perl -0pi -e 's/(#define\s+SMALLJAC_GENUS\s+)\d+/${1}2/' "$SMALLJAC_SRC/smalljac.h"
    grep -Eq '^[[:space:]]*#define[[:space:]]+SMALLJAC_GENUS[[:space:]]+2\b' "$SMALLJAC_SRC/smalljac.h" \
      || { echo "ERROR: failed to set SMALLJAC_GENUS=2" >&2; exit 1; }
  fi

  # --- build smalljac's lpdata against our ff_poly ---------------------------
  # Point INCLUDES/LIBS at our PREFIX instead of the makefile's default /usr/local.
  # Only the `lpdata` target is built (we don't need amicable/lpoly/moments).
  log "building smalljac lpdata ($SMALLJAC_TAR, genus 2)"
  make -C "$SMALLJAC_SRC" clean >/dev/null 2>&1 || true
  make -C "$SMALLJAC_SRC" -j"$JOBS" lpdata \
    CFLAGS="$SMALLJAC_CFLAGS" \
    INCLUDES="-I$PREFIX/include -I$PREFIX/include/ff_poly" \
    LIBS="-L$PREFIX/lib -lff_poly -lgmp -lm" >&2

  [ -x "$LPDATA_BIN" ] || { echo "ERROR: build finished but $LPDATA_BIN is missing" >&2; exit 1; }
fi

# --- smoke test: banner must report genus 2 (or compatible) -----------------
BANNER="$("$LPDATA_BIN" 2>&1 || true)"
if printf '%s' "$BANNER" | grep -qi 'genus'; then
  GENUS_LINE="$(printf '%s' "$BANNER" | grep -i 'genus' | head -1 | tr -s ' ')"
  log "lpdata banner: $GENUS_LINE"
fi

# --- optional: install into a shared bin (e.g. /usr/local/bin) --------------
if [ -n "${INSTALL_BIN:-}" ]; then
  log "installing lpdata into $INSTALL_BIN"
  mkdir -p "$INSTALL_BIN"
  cp -f "$LPDATA_BIN" "$INSTALL_BIN/lpdata"
  LPDATA_BIN="$INSTALL_BIN/lpdata"
fi

log "lpdata ready: $LPDATA_BIN"
# FINAL line of stdout = machine-parseable path (KEY=VALUE for easy `cut -d=`).
printf 'LPDATA=%s\n' "$LPDATA_BIN"
