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

static Lerror_t bad_supply(Lfunc *L)
{
  L->supply_ecode|=ERR_BAD_SUPPLY;
  return ERR_BAD_SUPPLY;
}

void Lfunc_use_lpoly(Lfunc_t Lf, uint64_t p, const acb_poly_t poly)
{
  Lfunc *L;
  L=(Lfunc *)Lf;
  if(L->raw_supplied) // can't push factors into raw-a_n overwrite mode
  {
    L->supply_ecode|=ERR_SUPPLY_CONFLICT; // void return: surfaced by Lfunc_compute
    return;
  }
  if(!poly)
  {
    bad_supply(L);
    return;
  }
  L->factor_supplied=true;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax(Lf);
    L->nmax_called=true;
  }
  use_lpoly(L,p,poly);
}

// Reduce the working coefficient count to new_M: we either ran out of data at,
// or were told to stop at, index/prime new_M+1. Also keep M0's direct block
// inside the reduced coefficient range, and keep buthe_M no larger than new_M.
// insufficient==true flags an unexpected shortfall (ERR_INSUFF_EULER); an
// explicit, trusted
// Lfunc_reduce_nmax passes false (a deliberate reduction is not an error).
// Returns the warning bit to OR into the caller's accumulator. This is the
// single funnel for both the callback zero-poly short-circuit and the array
// length-shortfall paths, so the M/M0/buthe_M clamp lives in exactly one place.
static Lerror_t shrink_M(Lfunc *L, uint64_t new_M, bool insufficient)
{
  L->M=new_M;
  if(new_M<UINT64_MAX && L->M0>new_M+1)
    L->M0=new_M+1; // direct block uses coefficients n < M0, so cap it at M
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
  if(L->raw_supplied) // raw a_n overwrote ans; the callback would multiply into it
  {
    L->supply_ecode|=ERR_SUPPLY_CONFLICT;
    return ERR_SUPPLY_CONFLICT;
  }
  L->factor_supplied=true;
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

// Shared body of the Euler-factor array front-ends: exactly one of fa (acb_poly
// array) and fz (fmpz_poly array) is non-NULL. We sieve the primes ourselves
// (exactly as Lfunc_use_all_lpolys) and hand the k-th factor to the shared
// use_lpoly for the k-th prime, so normalisation, inversion and the
// multiplicative sieve are not duplicated. An fmpz factor is converted into a
// single reused scratch acb_poly at working precision before use_lpoly, so
// Buthe's wf() still sees real per-prime forward and inverse factors. Running
// out of factors before nmax reduces M and warns, just like the callback
// zero-poly short-circuit; surplus factors (len > pi(nmax)) are ignored.
static Lerror_t use_lpolys_array(Lfunc *L, const acb_poly_struct *fa, const fmpz_poly_struct *fz, uint64_t len)
{
  if((fa && fz) || (!fa && !fz && len>0))
    return bad_supply(L);
  if(L->raw_supplied) // raw a_n overwrote ans; multiplying factors in is incoherent
  {
    L->supply_ecode|=ERR_SUPPLY_CONFLICT;
    return ERR_SUPPLY_CONFLICT;
  }
  L->factor_supplied=true;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax((Lfunc_t)L);
    L->nmax_called=true;
  }
  acb_poly_t g;
  acb_poly_init(g); // scratch for the fmpz->acb conversion (unused on the acb path)
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
    if(fa)
      use_lpoly(L,p,fa+k);
    else
    {
      acb_poly_set_fmpz_poly(g,fz+k,L->wprec); // converted at working precision
      use_lpoly(L,p,g);
    }
    k++;
  }
  primesieve_free_iterator(&it);
  acb_poly_clear(g);
  return ecode;
}

Lerror_t Lfunc_use_lpolys_acb(Lfunc_t Lf, const acb_poly_struct *f, uint64_t len)
{
  return use_lpolys_array((Lfunc *)Lf, f, NULL, len);
}

Lerror_t Lfunc_use_lpolys_fmpz(Lfunc_t Lf, const fmpz_poly_struct *f, uint64_t len)
{
  return use_lpolys_array((Lfunc *)Lf, NULL, f, len);
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

// Common fatal entry guards for the raw-a_n front-ends: (1) overwrite mode
// conflicts with any factor supply or a second raw supply; (2) the norm selector
// must be explicit and valid; (3) a_1 must be supplied and equal to 1. Also
// records the code in supply_ecode so Lfunc_compute bails even if the caller
// ignores the return value.
static Lerror_t raw_guard(Lfunc *L, uint64_t len, int norm_of_input, bool have_coeffs, bool a1_is_one)
{
  if(L->factor_supplied || L->raw_supplied)
  {
    L->supply_ecode|=ERR_SUPPLY_CONFLICT;
    return ERR_SUPPLY_CONFLICT;
  }
  if(norm_of_input!=ALGEBRAIC_NORM && norm_of_input!=ANALYTIC_NORM)
  {
    L->supply_ecode|=ERR_BAD_NORM;
    return ERR_BAD_NORM;
  }
  if(len==0)
  {
    L->supply_ecode|=ERR_A1_NOT_ONE;
    return ERR_A1_NOT_ONE;
  }
  if(!have_coeffs)
  {
    L->supply_ecode|=ERR_BAD_SUPPLY;
    return ERR_BAD_SUPPLY;
  }
  if(!a1_is_one)
  {
    L->supply_ecode|=ERR_A1_NOT_ONE;
    return ERR_A1_NOT_ONE;
  }
  return ERR_SUCCESS;
}

// Sanity-check a supplied analytic coefficient a_n against the degree's
// Euler-product Ramanujan bound |a_n| <= C*n^alpha (C, alpha set in g.c; the same
// bound the factor paths satisfy and M_error uses for the tail). Fatal
// ERR_COEFF_BOUND if |a_n| certainly exceeds it (arb_gt, so only a definite
// violation fires); catches a wrong normalisation_of_input or packed garbage.
static Lerror_t check_coeff_bound(Lfunc *L, const acb_t an, uint64_t n)
{
  int64_t prec=L->wprec;
  arb_t absn,bound,t;
  arb_init(absn);
  arb_init(bound);
  arb_init(t);
  acb_abs(absn,an,prec);
  arb_log_ui(t,n,prec);
  arb_mul(t,t,L->alpha,prec);
  arb_exp(t,t,prec);          // n^alpha
  arb_mul(bound,L->C,t,prec); // C * n^alpha
  Lerror_t e = arb_gt(absn,bound) ? ERR_COEFF_BOUND : ERR_SUCCESS;
  arb_clear(absn);
  arb_clear(bound);
  arb_clear(t);
  return e;
}

// Shared body of the raw Dirichlet-coefficient front-ends: exactly one of az
// (fmpz array) and aa (acb array) is non-NULL. The supplied a_n *overwrite*
// L->ans (the all-ones init), so they cannot be combined with any Euler-factor
// route. There are no per-prime factors, so RH verification is later skipped
// (gated on raw_supplied). The fmpz form writes each coefficient exactly; the
// acb form trusts the supplied ball. A short array reduces M and warns; surplus
// is ignored. The a_1 == 1 test (fmpz exact equality vs the acb ball containing
// 1) is the caller's, passed in as a1_is_one.
static Lerror_t use_dirichlet_coeffs(Lfunc *L, const fmpz *az, acb_srcptr aa, uint64_t len, int norm_of_input, bool have_coeffs, bool a1_is_one)
{
  Lerror_t guard=raw_guard(L,len,norm_of_input,have_coeffs,a1_is_one);
  if(guard)
    return guard;
  if(!L->nmax_called)
  {
    L->M=Lfunc_nmax((Lfunc_t)L);
    L->nmax_called=true;
  }
  L->raw_supplied=true;
  uint64_t use = (len<L->M) ? len : L->M;
  Lerror_t ecode=ERR_SUCCESS;
  for(uint64_t n=1;n<=use;n++)
  {
    if(az)
      acb_set_fmpz(L->ans[n-1],az+(n-1)); // exact
    else
      acb_set(L->ans[n-1],aa+(n-1)); // trust the supplied ball
    apply_input_norm(L->ans[n-1],n,norm_of_input,L);
    ecode|=check_coeff_bound(L,L->ans[n-1],n);
  }
  if(ecode)
    L->supply_ecode|=ecode; // surface a bound violation at compute time too
  if(len<L->M)
    ecode|=shrink_M(L,len,true);
  return ecode;
}

Lerror_t Lfunc_use_dirichlet_coeffs_fmpz(Lfunc_t Lf, const fmpz *a, uint64_t len, int norm_of_input)
{
  Lfunc *L=(Lfunc *)Lf;
  bool have_coeffs = len==0 || a!=NULL;
  bool a1_is_one = len>0 && a && fmpz_is_one(a+0);
  return use_dirichlet_coeffs(L,a,NULL,len,norm_of_input,have_coeffs,a1_is_one);
}

Lerror_t Lfunc_use_dirichlet_coeffs_acb(Lfunc_t Lf, acb_srcptr a, uint64_t len, int norm_of_input)
{
  Lfunc *L=(Lfunc *)Lf;
  // a_1's ball must contain 1 (certified contract: the ball encloses the truth)
  bool have_coeffs = len==0 || a!=NULL;
  bool a1_is_one = len>0 && a &&
    (arb_contains_si(acb_realref(a+0),1) && arb_contains_zero(acb_imagref(a+0)));
  return use_dirichlet_coeffs(L,NULL,a,len,norm_of_input,have_coeffs,a1_is_one);
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
