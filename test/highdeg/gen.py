#!/usr/bin/env python3
"""Generate the highdeg_check driver input for one object in objects.yaml.

Good Euler factors: lpdata (genus-1 EC a_p, then Sym^k) for kind ec/sympow;
Pari hyperellcharpoly (--backend pari) or genus-2 lpdata (--backend smalljac)
for kind genus2. Bad factors are injected verbatim from the YAML (LMFDB).
See test/highdeg/INTERFACES.md.

usage: gen.py LABEL [--objects FILE] [--driver PATH] [--lpdata PATH]
              [--backend pari|smalljac] [--workdir DIR]
"""
import argparse, atexit, math, os, shutil, subprocess, sys, tempfile
import yaml


def derive(obj):
    kind = obj["kind"]
    if kind == "ec":
        return dict(degree=2, norm=0.5, mus=[0, 1], sd=1)
    if kind == "genus2":
        return dict(degree=4, norm=0.5, mus=[0, 0, 1, 1], sd=1)
    if kind == "sympow":
        k = int(obj["sym"]); u = math.ceil(k / 2.0); mus = [0.0] * (k + 1)
        for i in range(u):
            mus[i] = -i; mus[i + u] = -i + 1
        if k % 2 == 0 and (k + 1) > k:
            mus[k] = -2 * math.floor(u / 2.0)
        return dict(degree=k + 1, norm=k / 2.0, mus=mus, sd=1)
    if kind == "cmf":
        # classical newform: degree-2*dim L-function (product of the dim Galois conjugates),
        # motivic weight k-1, dim Gamma_C(s+(k-1)/2) factors -> mus [0]*dim + [1]*dim.
        dim = int(obj["dim"])
        return dict(degree=2 * dim, norm=(int(obj["weight"]) - 1) / 2.0,
                    mus=[0] * dim + [1] * dim, sd=1)
    raise SystemExit("unknown kind: " + kind)


def query_nmax(driver, der, cond):
    with tempfile.NamedTemporaryFile("w", suffix=".hdr", delete=False) as f:
        f.write("%d %d %s %d\n%s\n" % (der["degree"], cond, der["norm"], der["sd"],
                                       " ".join(repr(x) for x in der["mus"])))
        hdr = f.name
    out = subprocess.run([driver, hdr], capture_output=True, text=True).stdout
    os.unlink(hdr)
    for line in out.splitlines():
        if line.startswith("nmax="):
            return int(line.split("=")[1])
    raise SystemExit("could not query nmax; driver said:\n" + out)


def lucas_sympoly(a, p, k):
    V = [2, a]
    for _ in range(2, k + 1):
        V.append(a * V[-1] - p * V[-2])
    poly = [1]
    def mul(P, Q):
        R = [0] * (len(P) + len(Q) - 1)
        for i, pi in enumerate(P):
            for j, qj in enumerate(Q):
                R[i + j] += pi * qj
        return R
    for j in range((k + 1) // 2):
        poly = mul(poly, [1, -(p ** j) * V[k - 2 * j], p ** k])
    if k % 2 == 0:
        poly = mul(poly, [1, -(p ** (k // 2))])
    return poly


def run_lpdata(lpdata, curvespec, nmax, workdir):
    pref = "lp"
    # flags=1 (SMALLJAC_GOOD_ONLY); NO jobs arg -> serial, single output file
    # "<pref>_lpdata.txt" (passing a jobs arg shards the filename in smalljac v4.1.3).
    subprocess.run([lpdata, pref, curvespec, str(nmax), "1"],
                   cwd=workdir, capture_output=True, text=True, check=True)
    pairs = []
    with open(os.path.join(workdir, pref + "_lpdata.txt")) as fh:
        for ln in fh:
            parts = ln.strip().split(",")
            if len(parts) < 2:
                continue
            try:
                pairs.append((int(parts[0]), [int(x) for x in parts[1:]]))
            except ValueError:
                continue  # header line
    return pairs


def good_ec_sympow(obj, der, nmax, lpdata, workdir, bad):
    k = int(obj.get("sym", 1))
    rows = []
    for p, fields in run_lpdata(lpdata, obj["curve"], nmax, workdir):
        if p > nmax or p in bad:
            continue
        a = -fields[0]  # lpdata t = -a_p
        rows.append((p, [1, -a, p] if k == 1 else lucas_sympoly(a, p, k)))
    return rows


def sextic_from_fh(f, h):
    h2 = [0] * (2 * len(h) - 1)
    for i, hi in enumerate(h):
        for j, hj in enumerate(h):
            h2[i + j] += hi * hj
    n = max(len(f), len(h2)); c = [0] * n
    for i, fi in enumerate(f): c[i] += 4 * fi
    for i, hi in enumerate(h2): c[i] += hi
    terms = []
    for i, ci in enumerate(c):
        if ci == 0:
            continue
        terms.append("%d" % ci if i == 0 else "%d*x" % ci if i == 1 else "%d*x^%d" % (ci, i))
    return "+".join(terms).replace("+-", "-") or "0"


def good_genus2(obj, der, nmax, backend, lpdata, workdir, bad):
    f, h = eval(obj["curve"])  # "[[f...],[h...]]"
    cond = obj["conductor"]
    rows = []
    if backend == "smalljac":
        for p, fields in run_lpdata(lpdata, sextic_from_fh(f, h), nmax, workdir):
            if p > nmax or p in bad or cond % p == 0:
                continue
            a1, a2 = fields[0], fields[1]
            rows.append((p, [1, a1, a2, a1 * p, p * p]))
        return rows
    # backend == pari: hyperellcharpoly via Sage
    script = (
        "from sage.all import primes, GF, PolynomialRing, pari\n"
        "f=%r; h=%r; cond=%d; nmax=%d\n"
        "for p in primes(nmax+1):\n"
        "    if cond %% p == 0: continue\n"
        "    Fp=GF(p); Rp=PolynomialRing(Fp,'x')\n"
        "    cp=pari.hyperellcharpoly([Rp(f),Rp(h)])\n"
        "    asc=[int(x) for x in cp.Vecrev()]\n"
        "    L=list(reversed(asc))\n"
        "    print(p, ' '.join(str(x) for x in L))\n" % (f, h, cond, nmax)
    )
    out = subprocess.run(["sage", "-python", "-c", script], capture_output=True, text=True, check=True).stdout
    for ln in out.splitlines():
        parts = ln.split()
        rows.append((int(parts[0]), [int(x) for x in parts[1:]]))
    return [(p, c) for (p, c) in rows if p not in bad and p <= nmax]  # injected factors win


def good_cmf(obj, der, nmax, bad):
    # Degree-2*dim factors of a classical newform's L-function (product over its dim Galois
    # conjugates). a_p come from Pari mfcoefs over the relative Hecke field K/Q(chi); the local
    # quadratic 1 - a_p X + eps(p) p^(k-1) X^2 is normed down to Q with two resultants (eliminate
    # the K/Q(chi) generator, then the character-field generator). eps(p) = (a_p^2 - a_{p^2})/p^(k-1)
    # is a character value keyed by p mod level. Good primes are p not dividing the level; the level's
    # prime divisors are the bad primes, supplied verbatim via bad_factors.
    N = int(obj["mf_level"]); k = int(obj["weight"])
    chi = int(obj["mf_chi"]); dim = int(obj["dim"])
    head = ("from sage.all import pari, euler_phi, prime_range\nimport sys\n"
            "N,k,chi,dim,nmax = %d,%d,%d,%d,%d\n" % (N, k, chi, dim, nmax))
    body = r'''
pari.default("parisizemax", "8G")                       # mfcoefs at large nmax needs headroom
cm = pari("Mod(%d,%d)" % (chi, N))
mf = pari.mfinit([N, k, cm], 0)
flds = pari.mffields(mf); B = pari.mfeigenbasis(mf)
cdeg = euler_phi(int(pari.znorder(cm)))                 # [Q(chi):Q]
idx = [i for i in range(len(flds)) if int(pari.poldegree(flds[i])) * cdeg == dim]
if len(idx) != 1:
    sys.exit("cmf: no unique eigenform of absolute degree %d (got %r)" % (dim, idx))
F = B[idx[0]]; Y = flds[idx[0]]; yv = Y.variable()
T = None                                                # character-field modulus
for j in range(int(pari.poldegree(Y)) + 1):
    c = pari.polcoef(Y, j, yv)
    if c.type() == "t_POLMOD":
        T = c.mod(); break
if T is None:
    sys.exit("cmf: trivial-character (rational) forms are not supported by this backend")
tv = T.variable(); X = pari("X")
co = pari.mfcoefs(F, nmax)
eps = {}
for p in prime_range(nmax + 1):
    if N % p == 0 or p * p > nmax:
        continue
    eps.setdefault(p % N, pari.lift((co[p] * co[p] - co[p * p]) / p**(k - 1)))
    if len(eps) == euler_phi(N):
        break
out = []
for p in prime_range(nmax + 1):
    if N % p == 0:
        continue
    loc = 1 - pari.lift(co[p]) * X + eps[p % N] * (p**(k - 1)) * X * X
    g = Y.polresultant(loc, yv)                         # norm K -> Q(chi): deg 2*dim/cdeg in X
    f8 = T.polresultant(pari.lift(g), tv)               # norm Q(chi) -> Q: deg 2*dim in X
    out.append("%d %s" % (p, " ".join(str(int(c)) for c in f8.Vecrev())))
sys.stdout.write("\n".join(out) + "\n")
'''
    out = subprocess.run(["sage", "-python", "-c", head + body],
                         capture_output=True, text=True, check=True).stdout
    rows = []
    for ln in out.splitlines():
        parts = ln.split()
        if parts:
            rows.append((int(parts[0]), [int(x) for x in parts[1:]]))
    return [(p, c) for (p, c) in rows if p not in bad and p <= nmax]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label")
    ap.add_argument("--objects", default=os.path.join(os.path.dirname(__file__), "objects.yaml"))
    ap.add_argument("--driver", default="build/test/highdeg_check.exe")
    ap.add_argument("--lpdata", default=os.environ.get("LPDATA", "/usr/local/bin/lpdata"))
    ap.add_argument("--backend", default=os.environ.get("BACKEND", "pari"))
    ap.add_argument("--workdir", default=None)
    a = ap.parse_args()
    # lpdata is invoked with cwd=workdir, so it MUST be an absolute path; abspath
    # the driver too so callers can pass relative paths.
    a.lpdata = os.path.abspath(a.lpdata)
    a.driver = os.path.abspath(a.driver)

    objs = yaml.safe_load(open(a.objects))["objects"]
    obj = next((o for o in objs if o["label"] == a.label), None)
    if obj is None:
        raise SystemExit("no object labelled " + a.label)
    der = derive(obj)
    cond = obj["conductor"]
    # `bad` = factors injected verbatim and skipped in good-factor generation:
    # real bad primes (bad_factors) PLUS force_factors (good primes the smalljac
    # model can't compute — e.g. p=2 for the completed-square genus-2 sextic).
    bad = {int(k): list(v) for k, v in (obj.get("bad_factors") or {}).items()}
    bad.update({int(k): list(v) for k, v in (obj.get("force_factors") or {}).items()})

    nmax = query_nmax(a.driver, der, cond)
    if a.workdir:
        workdir = a.workdir
    else:
        workdir = tempfile.mkdtemp(prefix="highdeg_")
        atexit.register(shutil.rmtree, workdir, ignore_errors=True)  # don't leak /tmp dirs

    if obj["kind"] in ("ec", "sympow"):
        rows = good_ec_sympow(obj, der, nmax, a.lpdata, workdir, bad)
    elif obj["kind"] == "cmf":
        rows = good_cmf(obj, der, nmax, bad)
    else:
        rows = good_genus2(obj, der, nmax, a.backend, a.lpdata, workdir, bad)

    for p, coeffs in bad.items():  # inject hardcoded bad factors
        if p <= nmax:
            rows.append((p, coeffs))

    e = obj["expected"]
    tay = e.get("taylor")
    expect = "EXPECT %d %s %s %s %s %s %s %d" % (
        e["rank"], repr(float(e["epsilon"][0])), repr(float(e["epsilon"][1])),
        str(e["first_zero"]), repr(float(e["first_zero_err"])),
        (str(tay) if tay is not None else "NA"),
        repr(float(e.get("taylor_err", 0.0))), 1 if e["tolerate_rh_error"] else 0)

    out = ["%d %d %s %d" % (der["degree"], cond, der["norm"], der["sd"]),
           " ".join(repr(x) for x in der["mus"]), expect]
    for p, coeffs in rows:
        coeffs = list(coeffs) + [0] * (der["degree"] + 1 - len(coeffs))
        out.append(str(p) + " " + " ".join(str(c) for c in coeffs))
    sys.stdout.write("\n".join(out) + "\n")
    sys.stderr.write("%s: kind=%s degree=%d cond=%d nmax=%d primes=%d backend=%s\n"
                     % (a.label, obj["kind"], der["degree"], cond, nmax, len(rows), a.backend))


if __name__ == "__main__":
    main()
