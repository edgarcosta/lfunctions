#include "inttypes.h"
#include "glfunc.h"
#include "glfunc_internals.h"
#include "primesieve.h"
#include <flint/acb_mat.h>

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

// Power / repeated-factor guard. Call at the top of Lfunc_compute, after the Euler
// factors have been supplied (so the signals are populated). Returns ERR_POWER if L
// looks like a perfect power / has a repeated primitive factor and the caller did not
// set allow_nonprimitive; ERR_SUCCESS otherwise.
//
// Belt-and-suspenders (sound, not complete): reject ONLY when BOTH
//   (1) no supplied full-degree local factor was proven squarefree, AND
//   (2) the empirical 2nd moment (1/#) sum |a_p|^2 >= POWER_MOMENT_THRESHOLD.
// (1) alone certifies a genuine primitive as safe, so a primitive is never false-rejected.
Lerror_t power_guard(Lfunc *L)
{
  if(L->allow_nonprimitive == YES)   // caller opted out of the guard
    return ERR_SUCCESS;
  if(L->seen_sqfree_fulldeg)         // rigorous certificate: no repeated factor
    return ERR_SUCCESS;
  if(L->moment_count == 0)           // no good primes seen; can't judge -> proceed
    return ERR_SUCCESS;
  arb_t S, thresh, diff;
  arb_init(S);
  arb_init(thresh);
  arb_init(diff);
  arb_div_ui(S, L->moment_sum, L->moment_count, L->wprec); // empirical 2nd moment
  arb_set_d(thresh, POWER_MOMENT_THRESHOLD);
  arb_sub(diff, S, thresh, L->wprec);
  bool power_like = arb_is_positive(diff); // S - threshold > 0 for the whole ball
  arb_clear(diff);
  arb_clear(thresh);
  arb_clear(S);
  return power_like ? ERR_POWER : ERR_SUCCESS;
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
