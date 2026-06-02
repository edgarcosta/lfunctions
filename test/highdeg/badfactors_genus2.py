#!/usr/bin/env python3
"""Compute a genus-2 curve's conductor and bad Euler factors via Pari lfungenus2.

REQUIRES A GIT/DEV BUILD OF PARI/GP (>= 2.18.1).  The RELEASED Pari (including
the version bundled with Sage, 2.17.1) has BUGGY genus-2 bad-Euler-factor code:
bug #2487, "lfungenus2 could give wrong result at odd primes" (it can silently
drop a bad local factor, e.g. returns [1] instead of [1,0,3]); fixed in 2.18.1.
MUST NOT be used here.  Build gp from https://pari.math.u-bordeaux.fr/git/pari.git
(./Configure --with-gmp && make gp) and point this script at it via --gp or $GP.

CAVEAT: Pari may omit the 2-adic part of the conductor (it prints N / 2^v with a
"conductor at 2" warning). For curves of even conductor, verify the prime-2 local
factor / conductor exponent separately. (All genus-2 curves currently in
objects.yaml have odd conductor, so this does not affect them.)

A curation helper for objects.yaml: smalljac (lpdata) supplies the GOOD factors,
the conductor + bad factors are hardcoded — this derives them straight from the
curve equation y^2 + h(x) y = f(x), so genus-2 objects need not already be in
LMFDB.  (For the LMFDB curves we use, the authoritative bad factors are taken
from g2c.bad_lfactors; this is the tool for curves outside LMFDB and for an
independent cross-check.)

usage:  GP=/path/to/git/gp test/highdeg/badfactors_genus2.py "[f0,f1,...]" "[h0,h1,...]"
        (ascending coefficients, as in objects.yaml's genus2 `curve: "[[f],[h]]"`)
"""
import argparse, os, subprocess, sys

ap = argparse.ArgumentParser()
ap.add_argument("f"); ap.add_argument("h")
ap.add_argument("--gp", default=os.environ.get("GP", "gp"))
a = ap.parse_args()
f, h = eval(a.f), eval(a.h)
spec = "[Polrev(%r),Polrev(%r)]" % (f, h)
gpcode = (
    "L=lfungenus2(%s);"
    "N=lfunparams(L)[1];"
    'print("conductor: ", N);'
    "fa=factor(N);"
    'for(i=1, matsize(fa)[1], p=fa[i,1]; print("  ", p, ": ", Vecrev(1/lfuneuler(L,p))));'
    "quit()"
) % spec
r = subprocess.run([a.gp, "-q"], input=gpcode, capture_output=True, text=True)
sys.stdout.write(r.stdout)
if r.returncode != 0:
    sys.stderr.write("gp failed (is --gp a git/dev build of Pari?):\n" + r.stderr)
    sys.exit(1)
