# High-degree regression suite — fixed interfaces (bead lfunctions-8dz)

These are the fixed contracts shared by `gen.py`, `highdeg_check.cpp`, and
`objects.yaml`; keep all three in sync when changing the format.

## A. `test/highdeg/objects.yaml` (data) — consumed by `gen.py`
```yaml
objects:
  - label: "11.a1"          # unique, filename-safe
    kind: ec                # ec | sympow | genus2
    curve: "[0,-1,1,-10,-20]"  # lpdata curve-spec: a-invariants for ec/sympow;
                               # for genus2: "[[f...],[h...]]" (y^2+h(x)y=f(x), ascending coeffs)
    sym: 1                  # symmetric power; REQUIRED for kind sympow; ignore for ec/genus2
    conductor: 11
    bad_factors:            # prime -> local L-poly coeffs ASCENDING, hardcoded from LMFDB
      11: [1, -1]
    expected:
      rank: 0
      epsilon: [1.0, 0.0]            # (re, im); for self-dual real objects this is +-1
      first_zero: "6.36261389471308870139"
      first_zero_err: 1.0e-15        # absolute half-width of the reference ball
      taylor: "0.2538418608..."      # OPTIONAL leading Taylor coeff (Lfunc_Taylor); real.
      taylor_err: 1.0e-12            #   include for ec/genus2 (LMFDB-golden, matches Lfunc_Taylor);
                                     #   OMIT for sympow (Pari/library leading-coeff normalisation
                                     #   differs) -> driver skips assertion 4 when taylor absent.
      tolerate_rh_error: false       # MUST be true for degree>=3 (see bead lfunctions-0zo)
```
`gen.py` DERIVES degree/normalisation/mus/self_dual from kind+sym (do NOT store them):
- ec     : degree 2, norm 0.5, mus [0,1], self_dual 1
- genus2 : degree 4, norm 0.5, mus [0,0,1,1], self_dual 1
- sympow : degree k+1, norm k/2, self_dual 1, mus by the ported formula:
    u=ceil(k/2); mus[i]=-i and mus[i+u]=-i+1 for i in [0,u); if k even and degree>k: mus[k]=-2*floor(u/2)

Optional object-level fields:
- `force_factors`: same shape as `bad_factors`; for GOOD primes the smalljac model cannot
  compute (p=2 of a completed-square genus-2 sextic, char 2). Injected verbatim and skipped in
  good-factor generation, exactly like `bad_factors`.
- `xfail: "reason"`: the suite RUNS the object but EXPECTS it to fail. A failure is reported
  `XFAIL` (suite stays green); an unexpected pass is `XPASS` (loud notice -> delete `xfail` to
  promote it to a real test). Used for the square-of-EC object (L = L(E)^2 has doubled zeros the
  solver cannot yet resolve).

## B. driver input file — emitted by `gen.py`, read by `test/highdeg_check.cpp`
```
<degree> <conductor> <normalisation> <self_dual:0|1>
<mu_0> <mu_1> ... <mu_{degree-1}>
EXPECT <rank> <eps_re> <eps_im> <z1> <z1_err> <taylor> <taylor_err> <tolerate_rh:0|1>
<p> <c0> <c1> ... <c_degree>      # local L-poly coeffs ASCENDING, decimal strings (bignum!)
...                                # one line per prime p <= Lfunc_nmax (good + bad), any order
```

## C. driver assertions (`highdeg_check.cpp`) — exit nonzero on ANY failure
1. `Lfunc_rank(L) == rank`
2. `acb_abs(Lfunc_sign(L))` overlaps 1  AND  `Lfunc_sign(L)` overlaps `(eps_re,eps_im)` widened by 1e-9
3. `Lfunc_zeros(L,0)[0]` overlaps `[z1 +- z1_err]`
4. `Lfunc_Taylor(L)` overlaps `[taylor +- taylor_err]`
5. fatal_error(ecode) must be false; if `tolerate_rh==0` also require `(ecode & ERR_RH_ERROR)==0`;
   if `tolerate_rh==1`, ignore ERR_RH_ERROR but no OTHER warning may be treated as failure.
Set `((Lfunc*)L)->self_dual = YES` when self_dual==1 (glfunc_internals.h). Coeffs via fmpz (int64 overflows for Sym^5+). Scrub `g_*` before each run.

## D. good-factor backends (in `gen.py`)
- ec / sympow: run `lpdata <prefix> "<curve>" <nmax> 1` (flag 1=good-only; no jobs arg — a jobs arg shards the output filename), parse `p,t` lines
  (t = linear coeff of EC L-poly = -a_p). For sympow apply Sym^k (Lucas V_m: V0=2,V1=a,Vm=a*V-p*V;
  factor = prod_{j=0..floor((k-1)/2)} (1 - p^j*V_{k-2j}*T + p^k*T^2), times (1 - p^{k/2}*T) if k even).
- genus2: `--backend pari` (LOCAL/validated): Pari `hyperellcharpoly([f_modp,h_modp])`, take `Vecrev`,
  reverse to get ascending L-poly. `--backend smalljac` (CI): genus-2 lpdata build (format TBD by the CI agent).
- bad primes: NEVER from lpdata; inject `bad_factors[p]` verbatim (pad to degree+1 with zeros).

Validated reference prototypes already in `bench/`: `sympow_bench.cpp` (driver core, has nmax-query +
self_dual + fmpz parsing), `gen.py` (sym powers), `gen_g2.py` (genus-2 via Pari, self-test OK).
