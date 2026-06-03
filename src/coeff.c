#include "inttypes.h"
#include "glfunc.h"
#include "glfunc_internals.h"
#include "primesieve.h"
#include <flint/acb_mat.h>
#include <flint/acb_poly.h>
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
        printf("Attempt to (re-)allocate memory for retained primes failed. Exiting.\n");
        exit(0);
      }
      acb_poly_struct *new_f = (acb_poly_struct *) realloc(L->retained_f, sizeof(acb_poly_struct)*nc);
      if (!new_f) {
        printf("Attempt to (re-)allocate memory for retained factors failed. Exiting.\n");
        exit(0);
      }
      L->retained_p = new_p;
      L->retained_f = new_f;
      L->retained_cap = nc;
    }
    L->retained_p[L->n_retained] = p;
    acb_poly_init(&L->retained_f[L->n_retained]);
    acb_poly_set(&L->retained_f[L->n_retained], f);
    L->n_retained++;
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

// root = f^(1/k) as an acb_poly (f must have unit constant term: good for L_p(0)=1).
static void poly_kth_root(acb_poly_t root, const acb_poly_t f, uint64_t k, int64_t prec)
{
  slong dlen = acb_poly_degree(f)/(slong)k + 2; // a couple of guard terms
  acb_t e; acb_init(e);
  acb_set_ui(e, 1); acb_div_ui(e, e, k, prec);   // e = 1/k
  acb_poly_pow_acb_series(root, f, e, dlen, prec);
  acb_clear(e);
}

// True iff f is rigorously a perfect k-th power of a degree deg(f)/k polynomial.
static bool poly_is_perfect_kth_power(const acb_poly_t f, uint64_t k, int64_t prec)
{
  slong d = acb_poly_degree(f);
  if (d < 1 || (d % (slong)k) != 0) return false;
  // The certificate is rigorous only for EXACT (point-ball) Euler factors: it proves
  // f = M_p^k by ball containment, meaningful only when f carries no width. All real
  // callers supply exact integer Euler factors; refuse to certify an inexact one so
  // we never claim a power we cannot prove.
  for (slong i = 0; i <= d; i++) {
    acb_srcptr ci = acb_poly_get_coeff_ptr(f, i);
    if (ci != NULL && !acb_is_exact(ci)) return false;
  }
  acb_poly_t root, chk; acb_poly_init(root); acb_poly_init(chk);
  poly_kth_root(root, f, k, prec);
  acb_poly_truncate(root, d/(slong)k + 1);    // drop spurious near-zero guard terms
  acb_poly_pow_ui(chk, root, (ulong)k, prec); // root^k, degree d
  // f being exact, root^k overlapping f certifies it: a non-power differs from any
  // k-th power by at least the integer coefficient gap (>= 1), far exceeding the
  // ~2^-prec width of the root^k ball, so an overlap can only mean f IS a k-th power.
  bool ok = acb_poly_overlaps(chk, f);
  acb_poly_clear(root); acb_poly_clear(chk);
  return ok;
}

// Fix and rigorously certify k for a pure power L = M^k. Returns ERR_SUCCESS with
// *k_out set when certified, ERR_POWER otherwise.
Lerror_t power_extract_prepare(Lfunc *L, uint64_t *k_out)
{
  if (L->moment_count == 0) return ERR_POWER;
  // candidate k from the 2nd moment (~ k^2)
  arb_t S; arb_init(S);
  arb_div_ui(S, L->moment_sum, L->moment_count, L->wprec);
  arb_sqrt(S, S, L->wprec);
  double kf = arf_get_d(arb_midref(S), ARF_RND_NEAR);
  arb_clear(S);
  uint64_t k = (uint64_t) (kf + 0.5);
  // A mis-read 2nd moment yields a safe false ERR_POWER: the rigorous gate certifies,
  // it cannot rescue a wrong candidate k.
  if (k < 2) return ERR_POWER;
  // necessary exact tell: conductor must be a perfect k-th power
  ulong rem, base = n_rootrem(&rem, (ulong)L->conductor, (ulong)k);
  (void)base;
  if (rem != 0) return ERR_POWER;
  // rigorous gate: every retained full-degree factor is a perfect k-th power
  bool seen_full = false;
  for (uint64_t i = 0; i < L->n_retained; i++) {
    if (acb_poly_degree(&L->retained_f[i]) == (slong)L->degree) {
      seen_full = true;
      if (!poly_is_perfect_kth_power(&L->retained_f[i], k, L->wprec))
        return ERR_POWER;
    }
  }
  if (!seen_full) return ERR_POWER;
  // The archimedean data must also be k copies: each sorted block of k mus must be
  // equal (mus are sorted in Lfunc_init_advanced). A mismatch means L's gamma factors
  // disagree with a pure power, so we must not extract (spec section 10).
  for (uint64_t i = 0; i + k <= L->degree; i += k)
    for (uint64_t j = 1; j < k; j++)
      if (L->mus[i] != L->mus[i + j])
        return ERR_POWER;
  *k_out = k;
  return ERR_SUCCESS;
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
  if(L->moment_count == 0)           // no good primes seen; can't judge -> proceed
    return ERR_SUCCESS;
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
