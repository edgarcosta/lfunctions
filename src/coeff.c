#include "inttypes.h"
#include "glfunc.h"
#include "glfunc_internals.h"
#include "primesieve.h"

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

void Lfunc_use_lpoly(Lfunc_t Lf, uint64_t p, const acb_poly_t poly)
{
  Lfunc *L;
  L=(Lfunc *)Lf;
  use_lpoly(L,p,poly);
}

// Reduce the working coefficient count to new_M: we either ran out of data at,
// or were told to stop at, index/prime new_M+1. Also keep buthe_M no larger than
// new_M. insufficient==true flags an unexpected shortfall (ERR_INSUFF_EULER); an
// explicit, trusted Lfunc_reduce_nmax passes false (a deliberate reduction is not
// an error). Returns the warning bit to OR into the caller's accumulator. This is
// the single funnel for both the callback zero-poly short-circuit and the array
// length-shortfall paths, so the M/buthe_M clamp lives in exactly one place.
static Lerror_t shrink_M(Lfunc *L, uint64_t new_M, bool insufficient)
{
  L->M=new_M;
#ifdef BUTHE
  if(L->buthe_M>new_M)
    L->buthe_M=new_M;
#endif
  return insufficient ? ERR_INSUFF_EULER : ERR_SUCCESS;
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
      ecode|=shrink_M(L,p-1,true); // we might get away with this
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

// Supply Euler factors as an array, one per consecutive prime (f[0] at p=2).
// We sieve the primes ourselves (exactly as Lfunc_use_all_lpolys) and hand the
// k-th factor to the shared use_lpoly for the k-th prime, so normalisation,
// inversion and the multiplicative sieve are not duplicated. Running out of
// factors before nmax reduces M and warns, just like the callback zero-poly
// short-circuit; surplus factors (len > pi(nmax)) are ignored.
Lerror_t Lfunc_use_lpolys_acb(Lfunc_t Lf, const acb_poly_struct *f, uint64_t len)
{
  Lfunc *L=(Lfunc *)Lf;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t p=0, k=0;
  Lerror_t ecode=ERR_SUCCESS;
  while((p=primesieve_next_prime(&it)) <= L->M)
  {
    if(k>=len) // ran out of supplied factors before nmax
    {
      ecode|=shrink_M(L,p-1,true);
      break;
    }
    use_lpoly(L,p,f+k);
    k++;
  }
  primesieve_free_iterator(&it);
  return ecode;
}

// As Lfunc_use_lpolys_acb, but the factors are exact integer polynomials. Each
// is converted to a reused scratch acb_poly (exact for coefficients that fit in
// the working precision, which Ramanujan-bounded local factors always do) and
// passed to the identical use_lpoly path, so Buthe's wf() still sees real
// per-prime forward and inverse factors.
Lerror_t Lfunc_use_lpolys_fmpz(Lfunc_t Lf, const fmpz_poly_struct *f, uint64_t len)
{
  Lfunc *L=(Lfunc *)Lf;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }
  acb_poly_t g;
  acb_poly_init(g);
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t p=0, k=0;
  Lerror_t ecode=ERR_SUCCESS;
  while((p=primesieve_next_prime(&it)) <= L->M)
  {
    if(k>=len) // ran out of supplied factors before nmax
    {
      ecode|=shrink_M(L,p-1,true);
      break;
    }
    acb_poly_set_fmpz_poly(g,f+k,L->wprec);
    use_lpoly(L,p,g);
    k++;
  }
  primesieve_free_iterator(&it);
  acb_poly_clear(g);
  return ecode;
}

// Declare the Ramanujan-type growth bound |a_n| <= C * n^alpha on the *analytic*
// coefficients. The raw-a_n paths feed these to M_error in place of the
// Euler-product default that g.c installs (alpha=1, C=coeff_bound(degree)); only
// the caller can certify the unsupplied tail beyond the supplied coefficients.
// compute_g runs at init (before this), so the values stored here survive to
// Lfunc_compute; the g.c assignment is also guarded on coeff_bound_set as a
// belt-and-suspenders against that ordering ever changing.
Lerror_t Lfunc_set_coeff_bound(Lfunc_t Lf, const arb_t C, double alpha)
{
  Lfunc *L=(Lfunc *)Lf;
  arb_set(L->C,C);            // C/alpha were arb_init'd by compute_g at init
  arb_set_d(L->alpha,alpha);
  L->coeff_bound_set=true;
  return ERR_SUCCESS;
}

// Move one supplied coefficient at index n (1-based) into the analytic
// normalisation: ALGEBRAIC_NORM multiplies by n^{-normalisation} (the direct
// analogue of use_lpoly's per-factor p^{-m*normalisation}); ANALYTIC_NORM, the
// n=1 term, and normalisation 0 need no shift.
static void apply_input_norm(acb_t z, uint64_t n, int norm_of_input, Lfunc *L)
{
  if(norm_of_input==ANALYTIC_NORM || n==1 || L->normalisation==0.0)
    return;
  int64_t prec=L->wprec;
  arb_t logn,f;
  arb_init(logn);
  arb_init(f);
  arb_log_ui(logn,n,prec);
  arb_set_d(f,-L->normalisation);
  arb_mul(f,f,logn,prec);
  arb_exp(f,f,prec);          // n^{-normalisation}
  acb_mul_arb(z,z,f,prec);
  arb_clear(logn);
  arb_clear(f);
}

// Supply the Dirichlet coefficients a_n directly (a[0]=a_1). These *overwrite*
// L->ans (the all-ones init), so they cannot be combined with any Euler-factor
// route. There are no per-prime factors, so RH verification is later skipped
// (no_lpolys). A short array reduces M and warns; surplus is ignored.
Lerror_t Lfunc_use_dirichlet_coeffs_fmpz(Lfunc_t Lf, const fmpz *a, uint64_t len, int norm_of_input)
{
  Lfunc *L=(Lfunc *)Lf;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }
  L->raw_supplied=true;
  L->no_lpolys=true;
  uint64_t use = (len<L->M) ? len : L->M;
  for(uint64_t n=1;n<=use;n++)
  {
    acb_set_fmpz(L->ans[n-1],a+(n-1)); // exact
    apply_input_norm(L->ans[n-1],n,norm_of_input,L);
  }
  if(len<L->M)
    return shrink_M(L,len,true);
  return ERR_SUCCESS;
}

Lerror_t Lfunc_use_dirichlet_coeffs_acb(Lfunc_t Lf, acb_srcptr a, uint64_t len, int norm_of_input)
{
  Lfunc *L=(Lfunc *)Lf;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }
  L->raw_supplied=true;
  L->no_lpolys=true;
  uint64_t use = (len<L->M) ? len : L->M;
  for(uint64_t n=1;n<=use;n++)
  {
    acb_set(L->ans[n-1],a+(n-1)); // trust the supplied ball
    apply_input_norm(L->ans[n-1],n,norm_of_input,L);
  }
  if(len<L->M)
    return shrink_M(L,len,true);
  return ERR_SUCCESS;
}

bool Lfunc_reduce_nmax(Lfunc_t LL, uint64_t nmax)
{
  Lfunc *L=(Lfunc *)LL;
  uint64_t M;
  M=Lfunc_nmax(LL); // what is the current M
  if(nmax>=M) // I won't let you increase it
    return false;
  shrink_M(L,nmax,false); // trusted, explicit reduction: not a shortfall
  return true;
}


#ifdef __cplusplus
}
#endif
