#include "inttypes.h"
#include "glfunc.h"
#include "glfunc_internals.h"
#include "primesieve.h"
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>

#ifdef __cplusplus
extern "C"{
#endif



uint64_t Lfunc_nmax(Lfunc_t Lf)
{
  Lfunc *L;
  L=(Lfunc *)Lf;

  if(L->nmax_called)
    return L->M;

  int64_t prec=L->wprec;
  arb_t tmp,tmp1;
  arb_init(tmp);arb_init(tmp1);
  arb_sqrt_ui(tmp,L->conductor,prec);
  arb_inv(L->one_over_root_N,tmp,prec);
  complete_ftwiddle_error(L,prec);
  if(verbose){printf("Final Ftwiddle Error set to ");arb_printd(L->ftwiddle_error,10);printf("\n");}

  L->dc=sqrt((double) L->conductor);
  L->M0=ceil(L->dc/100);
  if(verbose)printf("M0 set to %" PRIu64 ".\n",L->M0);
  L->M=L->dc*exp(2*M_PI*(L->hi_i+0.5)*L->one_over_B);
  if(verbose)printf("M computed from hi_i = %" PRIu64 "\n",L->M);

  /*
  // I think this attempt to reduce M empirically is not worth the effort
  arb_zero(tmp1);
  uint64_t old_M=L->M;
  while(true)
  {
  L->M=(double)L->M/1.05;
  M_error(tmp,tmp1,L,prec);
  arb_mul_2exp_si(tmp,tmp,L->target_prec+75);
  arb_sub_ui(tmp,tmp,1,prec);
  if(!arb_is_negative(tmp))
  break;
  old_M=L->M;
  }
  L->M=old_M;
  if(verbose)printf("M reduced to %" PRIu64 ".\n",L->M);
  */

  if(L->M>L->allocated_M)
  {
    //printf("Need more space for Dirichlet coefficients.\n");
    for(size_t i = 0; i < L->allocated_M; ++i)
      acb_clear(L->ans[i]);

    while(L->allocated_M < L->M) {
      L->allocated_M<<=1;
      if(L->allocated_M == 0) // L->M was huge!
        L->allocated_M = L->M;
    }
    L->ans=(acb_t *)realloc(L->ans,sizeof(acb_t)*L->allocated_M);
    if(!L->ans)
    {
      printf("Attempt to (re-)allocate memory for Dirichlet coefficients failed. Exiting.\n");
      exit(0);
    }
    for(size_t i = 0; i < L->allocated_M; ++i)
      acb_init(L->ans[i]);

    //printf("re-allocated enough memory for Dirichlet coefficients.\n");
  }

  for(size_t i = 0; i < L->M; ++i)
    acb_set_ui(L->ans[i],1);
#ifdef BUTHE
  arb_zero(L->buthe_Wf);
  L->buthe_M=sqrt((double) L->M);
#endif
  L->nmax_called=true;

  arb_clear(tmp);arb_clear(tmp1);
  return L->M;
}

// Rigorously decide whether the polynomial f (degree d >= 1) is squarefree
// (distinct roots), via res(f, f') computed as the Sylvester-matrix determinant.
// Returns true ONLY if the determinant ball provably excludes zero, so a `true`
// is a rigorous certificate that f has no repeated root. A `false` means either
// f has a repeated root or the balls were too wide to decide -- both are safe to
// treat as "not certified squarefree".
static bool poly_is_squarefree_certified(const acb_poly_t f, int64_t prec)
{
  slong d = acb_poly_degree(f);
  if(d < 1)
    return false;
  acb_poly_t fp;
  acb_poly_init(fp);
  acb_poly_derivative(fp, f, prec);     // f', degree d-1
  slong n = 2*d - 1;                    // Sylvester matrix of (f, f') is n x n
  acb_mat_t S;
  acb_mat_init(S, n, n);
  acb_mat_zero(S);
  // top d-1 rows: coeffs of f (highest degree first), row i shifted right by i
  for(slong i = 0; i < d-1; i++)
    for(slong j = 0; j <= d; j++)
      acb_poly_get_coeff_acb(acb_mat_entry(S, i, i+j), f, d-j);
  // bottom d rows: coeffs of f' (highest degree first), row i shifted right by i
  for(slong i = 0; i < d; i++)
    for(slong j = 0; j <= d-1; j++)
      acb_poly_get_coeff_acb(acb_mat_entry(S, d-1+i, i+j), fp, (d-1)-j);
  acb_t det;
  acb_init(det);
  acb_mat_det(det, S, prec);
  bool sqfree = !acb_contains_zero(det);
  acb_clear(det);
  acb_mat_clear(S);
  acb_poly_clear(fp);
  return sqfree;
}

#ifdef BUTHE
void use_inv_lpoly(Lfunc *L, uint64_t p, acb_poly_t c, acb_poly_t f, uint64_t prec)
  #else
void use_inv_lpoly(Lfunc *L, uint64_t p, acb_poly_t c, uint64_t prec)
#endif  
{
  acb_t tmp;
  acb_init(tmp);
  //if(p==2) {printf("p=%" PRIu64 " 1/poly=",p);acb_poly_printd(c,20);printf("\npoly=");acb_poly_printd(f,20);printf("\n");}
  #ifdef BUTHE
  wf(L, p, c, f, prec); // do the Buthe bit, see buthe.c
  #endif
  // use inverted poly to populate Dirichlet coefficients
  uint64_t pnn=p*p, pn=p,pow=1;
  while(pn <= L->M) {
    acb_poly_get_coeff_acb(tmp, c, pow);
    uint64_t ptr = pn, count = 1;
    while(ptr <= L->M) {
      if(count < p) {// its not a higher prime power
        acb_mul(L->ans[ptr-1], L->ans[ptr-1], tmp, prec);
        count++;
        ptr += pn;
      } else {// it is higher prime power, so skip it
        ptr += pn;
        count = 1;
      }
    }
    pn *= p;
    pnn *= p;
    pow++;
  }
  acb_clear(tmp);
}

void use_lpoly(Lfunc *L, uint64_t p, const acb_poly_t f)
{
  int64_t prec=L->wprec;
  acb_t tmp;
  arb_t logp,tmp1,tmp2;
  acb_poly_t n_poly,inv_poly;
  arb_init(logp);
  acb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  acb_poly_init(n_poly);
  acb_poly_init(inv_poly);
  // w70.2: retain the raw Euler factor for later k-th-root extraction
  if (L->extract_powers == YES) {
    if (L->n_retained == L->retained_cap) {
      uint64_t nc = L->retained_cap ? 2*L->retained_cap : 256;
      uint64_t *new_p = (uint64_t *) realloc(L->retained_p, sizeof(uint64_t)*nc);
      if (!new_p) {
        L->extract_powers = NO;
      } else {
        acb_poly_struct *new_f = (acb_poly_struct *) realloc(L->retained_f, sizeof(acb_poly_struct)*nc);
        if (!new_f) {
          // Note: new_p doesn't strictly leak here since we realloc'd in-place or returned new pointer,
          // but we assign it to retained_p so it gets freed on clear.
          L->retained_p = new_p; 
          L->extract_powers = NO;
        } else {
          L->retained_p = new_p;
          L->retained_f = new_f;
          L->retained_cap = nc;
        }
      }
    }
    if (L->extract_powers == YES) {
      L->retained_p[L->n_retained] = p;
      acb_poly_init(&L->retained_f[L->n_retained]);
      acb_poly_set(&L->retained_f[L->n_retained], f);
      L->n_retained++;
    }
  }
  //if(p<=2){printf("in use_lpoly pre-norm with p = %" PRIu64 "\n",p);acb_poly_printd(f,20);printf("\n");}
  arb_log_ui(logp,p,prec);
  // normalise by multiplying each term by p^(-m norm)
  acb_poly_one(n_poly);
  for(int64_t m=1;m<acb_poly_length(f);m++)
  {
    arb_set_d(tmp1,-L->normalisation);
    arb_mul(tmp2,tmp1,logp,prec);
    arb_mul_ui(tmp1,tmp2,m,prec);
    arb_exp(tmp2,tmp1,prec);
    acb_mul_arb(tmp,acb_poly_get_coeff_ptr(f,m),tmp2,prec);
    acb_poly_set_coeff_acb(n_poly,m,tmp);
  }
  //if(p<=11){printf("in use_lpoly post-norm with p = %" PRIu64 "\n",p);acb_poly_printd(n_poly,20);printf("\n");}
  uint64_t k=1,pk=p;
  while(pk<=L->M) {k++;pk*=p;}

  acb_poly_inv_series(inv_poly,n_poly,k,prec);
  //if(p<=11){printf("Inverted poly\n");
  //acb_poly_printd(inv_poly,20);printf("\n------------------\n");
  //}
  // ---- power-guard detection signals, accumulated over full-degree (good) primes ----
  // a_p (analytic) is the T^1 coefficient of the inverse local factor.
  if(acb_poly_degree(f) == (slong)L->degree)
  {
    acb_t ap;
    arb_t aap;
    acb_init(ap);
    arb_init(aap);
    acb_poly_get_coeff_acb(ap, inv_poly, 1);
    acb_abs(aap, ap, prec);            // |a_p|
    arb_mul(aap, aap, aap, prec);      // |a_p|^2
    arb_add(L->moment_sum, L->moment_sum, aap, prec);
    L->moment_count++;
    // one squarefree full-degree factor rigorously certifies "no repeated factor"
    if(!L->seen_sqfree_fulldeg && L->moment_count <= POWER_SQFREE_PROBES)
      if(poly_is_squarefree_certified(f, prec))
        L->seen_sqfree_fulldeg = true;
    arb_clear(aap);
    acb_clear(ap);
  }
  #ifdef BUTHE
  use_inv_lpoly(L,p,inv_poly,n_poly,prec);
#else
  use_inv_lpoly(L,p,inv_poly,prec);
#endif
  arb_clear(logp);
  acb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  acb_poly_clear(n_poly);
  acb_poly_clear(inv_poly);

}

// Cheap, exact "is L a pure power?" tell from the conductor alone: cond(M^k) = cond(M)^k,
// so a pure power has a perfect-k-th-power conductor. Necessary, not sufficient (a primitive
// can have a perfect-power conductor), so this only ever corroborates -- Signal 1 (the
// squarefree certificate) is what keeps the guard from false-rejecting such primitives.
static bool conductor_is_perfect_power(uint64_t N)
{
  if(N < 4)            // 0,1,2,3 are not (nontrivial) perfect powers
    return false;
  ulong root;
  return n_is_perfect_power(&root, (ulong) N) != 0;
}

// Exact k-th root of the conductor: returns true iff `cond` is a perfect k-th power
// (cond = base^k for some integer base), setting *base to that root. This is the
// "is it THIS power" question (k fixed), distinct from conductor_is_perfect_power's
// "is it ANY power". Used by both power_extract_prepare and extract_and_assemble.
bool conductor_kth_root(uint64_t cond, uint64_t k, uint64_t *base)
{
  ulong rem, b = n_rootrem(&rem, (ulong)cond, (ulong)k);
  *base = (uint64_t)b;
  return rem == 0;
}

// Rigorously decide whether the EXACT integer Euler factor Lp is a perfect k-th
// power of an algebraic polynomial. If so, set Mp to that algebraic root and return
// true; else return false. EXACT: rounds the ball k-th root to algebraic integers and checks
// Mp^k == Lp exactly over acb.
bool poly_exact_kth_root(acb_poly_t Mp, const acb_poly_t Lp, uint64_t k, int64_t prec)
{
  slong d = acb_poly_degree(Lp);
  if (d < 0 || (d % (slong)k) != 0) return false;

  // We require Lp to be exact
  for (slong i = 0; i <= d; i++) {
    acb_srcptr c = acb_poly_get_coeff_ptr(Lp, i);
    if (c != NULL && !acb_is_exact(c)) return false;
  }

  acb_poly_t root; acb_poly_init(root);
  acb_t e; acb_init(e); acb_set_ui(e, 1); acb_div_ui(e, e, k, prec);
  acb_poly_pow_acb_series(root, Lp, e, d/(slong)k + 1, prec);
  acb_clear(e);

  acb_poly_zero(Mp);
  for (slong i = 0; i <= d/(slong)k; i++) {
    acb_srcptr c = acb_poly_get_coeff_ptr(root, i);
    acb_t z; acb_init(z);
    if (c != NULL) {
      fmpz_t r, im; fmpz_init(r); fmpz_init(im);
      arf_get_fmpz(r, arb_midref(acb_realref(c)), ARF_RND_NEAR);
      arf_get_fmpz(im, arb_midref(acb_imagref(c)), ARF_RND_NEAR);
      arb_set_fmpz(acb_realref(z), r);
      arb_set_fmpz(acb_imagref(z), im);
      fmpz_clear(r); fmpz_clear(im);
    }
    acb_poly_set_coeff_acb(Mp, i, z);
    acb_clear(z);
  }
  acb_poly_clear(root);

  acb_poly_t chk; acb_poly_init(chk);
  acb_poly_pow_ui(chk, Mp, (ulong)k, prec);

  bool ok = true;
  for (slong i = 0; i <= d && ok; i++) {
    acb_srcptr c_L = acb_poly_get_coeff_ptr(Lp, i);
    acb_srcptr c_chk = acb_poly_get_coeff_ptr(chk, i);
    acb_t exact_chk; acb_init(exact_chk);
    if (c_chk != NULL) {
      fmpz_t r, im; fmpz_init(r); fmpz_init(im);
      arf_get_fmpz(r, arb_midref(acb_realref(c_chk)), ARF_RND_NEAR);
      arf_get_fmpz(im, arb_midref(acb_imagref(c_chk)), ARF_RND_NEAR);
      arb_set_fmpz(acb_realref(exact_chk), r);
      arb_set_fmpz(acb_imagref(exact_chk), im);
      fmpz_clear(r); fmpz_clear(im);
    }
    if (c_L == NULL) {
      if (!acb_is_zero(exact_chk)) ok = false;
    } else {
      if (!acb_equal(c_L, exact_chk)) ok = false;
    }
    acb_clear(exact_chk);
  }
  acb_poly_clear(chk);
  return ok;
}

// Free and reset the certified-root store (the partial work of a rejected candidate k,
// or any leftover from a prior call). After this the store is empty and owns nothing.
static void cert_roots_reset(Lfunc *L)
{
  if (L->cert_roots) {
    for (uint64_t i = 0; i < L->n_cert_roots; i++)
      acb_poly_clear(&L->cert_roots[i]);
    free(L->cert_roots);
    L->cert_roots = NULL;
  }
  L->n_cert_roots = 0;
}

// Run the full rigorous certificate for a fixed candidate k: EVERY retained factor
// (good AND bad) must be an exact integer k-th power (certified over acb_poly), each
// sorted block of k mus must be equal, and at least one full-degree (good) factor must
// be present. On success, the exact algebraic roots are kept in L->cert_roots (one per
// retained factor, same order) for reuse when feeding M, and true is returned. On any
// failure the partially-built store is freed and false is returned.
static bool power_certify_k(Lfunc *L, uint64_t k)
{
  // The archimedean data must be k copies: each sorted block of k mus must be equal
  // (mus are sorted in Lfunc_init_advanced). A mismatch means L's gamma factors disagree
  // with a pure power, so we must not extract (spec section 10). Cheap; check first.
  for (uint64_t i = 0; i + k <= L->degree; i += k)
    for (uint64_t j = 1; j < k; j++)
      if (L->mus[i] != L->mus[i + j])
        return false;
  // rigorous per-prime gate: a single uncertified factor would make L = M^k false, so the
  // k-th root fed to M would be unsound; reject. Keep each certified root for reuse.
  cert_roots_reset(L);
  if (L->n_retained) {
    L->cert_roots = (acb_poly_struct *) malloc(sizeof(acb_poly_struct)*L->n_retained);
    if (!L->cert_roots) {
      return false;
    }
  }
  bool seen_full = false;
  for (uint64_t i = 0; i < L->n_retained; i++) {
    if (acb_poly_degree(&L->retained_f[i]) == (slong)L->degree) seen_full = true;
    acb_poly_init(&L->cert_roots[i]);
    L->n_cert_roots = i + 1;
    if (!poly_exact_kth_root(&L->cert_roots[i], &L->retained_f[i], k, L->wprec)) {
      cert_roots_reset(L);
      return false;
    }
  }
  if (!seen_full) { cert_roots_reset(L); return false; }
  return true;
}

// Fix and rigorously certify k for a pure power L = M^k. Returns ERR_SUCCESS with
// *k_out set when certified, ERR_POWER otherwise. On success the certified integer
// roots are cached in L->cert_roots for reuse by extract_and_assemble.
Lerror_t power_extract_prepare(Lfunc *L, uint64_t *k_out)
{
  if (L->moment_count == 0) return ERR_POWER;
  
  // For L = M^k the true exponent k satisfies k | degree (M has degree/k Gamma factors)
  // and makes the conductor an exact k-th power (cond(M^k) = cond(M)^k). 
  // We consider every k in [2, degree] with k | degree and cond an exact k-th power,
  // preferring the LARGEST that passes the rigorous certificate. 
  // (We test conductor-exactness per candidate via conductor_kth_root rather
  // than n_is_perfect_power, whose returned exponent is not necessarily maximal -- it gives
  // 2 for 37^4 = 1369^2.) The per-prime exact gate and the mus-block gate still run
  // for the chosen k and reject anything unsound, so a too-large k is filtered back down
  // (e.g. cond a 4th power but the factors only squares).
  uint64_t cbase;
  for (uint64_t k = L->degree; k >= 2; k--) {
    if ((L->degree % k) != 0) continue;   // M's degree = degree/k must be a positive integer
    if (!conductor_kth_root(L->conductor, k, &cbase)) continue; // cond not an exact k-th power
    if (power_certify_k(L, k)) { *k_out = k; return ERR_SUCCESS; }
  }
  cert_roots_reset(L);
  return ERR_POWER;
}

// Power / repeated-factor guard. Call at the top of Lfunc_compute, after the Euler
// factors have been supplied (so the signals are populated). Returns ERR_POWER if L
// looks like a perfect power / has a repeated primitive factor, ERR_SUCCESS if it
// looks primitive. The verdict is independent of extract_powers; Lfunc_compute then
// rejects (opted out) or extracts-and-assembles (opted in) on an ERR_POWER verdict.
//
// Belt-and-suspenders (sound, not complete): reject ONLY when (1) is false, AND (2a OR 2b):
//   (1)  some supplied full-degree local factor was proven squarefree, OR
//   (2a) the empirical 2nd moment (1/#) sum |a_p|^2 >= POWER_MOMENT_THRESHOLD, OR
//   (2b) the conductor is a perfect power (cond(M^k)=cond(M)^k, exact tell for PURE powers).
// (1) alone certifies a genuine primitive as safe and is checked FIRST, so a primitive is
// never false-rejected even if it happens to have a perfect-power conductor.
Lerror_t power_guard(Lfunc *L)
{
  // Detection is independent of extract_powers: Lfunc_compute consults this verdict
  // and then either rejects (opted out) or extracts-and-assembles (opted in).
  if(L->seen_sqfree_fulldeg)         // rigorous certificate: no repeated factor
    return ERR_SUCCESS;
  if(L->moment_count == 0)           // no good primes seen; can't judge from moments.
    // A perfect-power conductor with NO good primes is the degenerate case we cannot
    // certify primitive (Signal 1 didn't fire either): flag it rather than silently
    // proceed (which could feed doubled zeros into the pipeline). A genuine primitive
    // with a squarefree good prime is already handled by the check above.
    return conductor_is_perfect_power(L->conductor) ? ERR_POWER : ERR_SUCCESS;
  // Signal 2a (heuristic): empirical 2nd moment >= threshold.
  arb_t S, thresh, diff;
  arb_init(S);
  arb_init(thresh);
  arb_init(diff);
  arb_div_ui(S, L->moment_sum, L->moment_count, L->wprec);
  arb_set_d(thresh, POWER_MOMENT_THRESHOLD);
  arb_sub(diff, S, thresh, L->wprec);
  bool moment_like = arb_is_positive(diff); // S - threshold > 0 for the whole ball
  arb_clear(diff);
  arb_clear(thresh);
  arb_clear(S);
  // Signal 2b (cheap exact tell for PURE powers): cond is a perfect k-th power.
  bool conductor_like = conductor_is_perfect_power(L->conductor);
  return (moment_like || conductor_like) ? ERR_POWER : ERR_SUCCESS;
}

void Lfunc_use_lpoly(Lfunc_t Lf, uint64_t p, const acb_poly_t poly)
{
  Lfunc *L;
  L=(Lfunc *)Lf;
  use_lpoly(L,p,poly);
}



// call lpoly_callback with every prime <=L->M
// so we can populate ans, the Dirichlet coefficients
Lerror_t Lfunc_use_all_lpolys(Lfunc_t Lf, void (*lpoly_callback) (acb_poly_t lpoly, uint64_t p, int d, int64_t prec, void *parm), void *param)
{
  Lfunc *L;
  L=(Lfunc *)Lf;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }

  acb_poly_t lp;
  acb_poly_init(lp);
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t p=0;
  Lerror_t ecode=ERR_SUCCESS;
  while((p=primesieve_next_prime(&it)) <= L->M)
  {
    lpoly_callback(lp,p,L->degree,L->wprec,param);
    if(acb_poly_is_zero(lp)) // ran out of Euler polys
    {
      #ifdef BUTHE
      if(p<L->buthe_M)
        L->buthe_M=p-1; // this is likely to mean we compute garbage
      #endif
      L->M=p-1; // we might get away with this
      ecode|=ERR_INSUFF_EULER;
      break;
    }
    use_lpoly(L,p,lp);
  }

  //for(i=0;i<20;i++)
  //{printf("Coefficient %" PRIu64 " set to ",i+1);acb_printd(L->ans[i],20);printf("\n");}

  primesieve_free_iterator(&it);
  acb_poly_clear(lp);
  return ecode;
}

bool Lfunc_reduce_nmax(Lfunc_t LL, uint64_t nmax)
{
  Lfunc *L=(Lfunc *)LL;
  uint64_t M;
  M=Lfunc_nmax(LL); // what is the current M
  if(nmax>=M) // I won't let you increase it
    return false;
  L->M=nmax;
  #ifdef BUTHE
  if(L->buthe_M>nmax) // we could be in serious trouble here
    L->buthe_M=nmax;
  #endif
  return true;
}


#ifdef __cplusplus
}
#endif
