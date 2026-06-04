#include <flint/acb_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/ulong_extras.h>
#include "glfunc.h"
#include "glfunc_internals.h"
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

// only used for debugging purposes
void acb_print_max_abs_error(acb_t *vec, uint64_t n) {
  arb_t tmp,tmp1;
  arb_init(tmp);arb_init(tmp1);

  arb_get_rad_arb(tmp,acb_realref(vec[0]));
  arb_get_rad_arb(tmp1,acb_imagref(vec[0]));
  if(arb_gt(tmp1,tmp))
    arb_set(tmp,tmp1);
  for(uint64_t nn=1;nn<n;nn++)
  {
    arb_get_rad_arb(tmp1,acb_realref(vec[nn]));
    if(arb_gt(tmp1,tmp))
      arb_set(tmp,tmp1);
    arb_get_rad_arb(tmp1,acb_imagref(vec[nn]));
    if(arb_gt(tmp1,tmp))
      arb_set(tmp,tmp1);
  }
  arb_printd(tmp,20);
  arb_clear(tmp);
  arb_clear(tmp1);
}

void acb_print_max_rel_error(acb_t *vec, uint64_t n)
{
  int64_t wc=acb_rel_error_bits(vec[0]);
  for(uint64_t i=1;i<n;i++)
  {
    int64_t e=acb_rel_error_bits(vec[i]);
    //acb_printd(vec[i],20);printf(" %" PRId64 "\n",e);
    if(e>wc)
      wc=e;
  }
  printf("%" PRId64,wc);
}

// copy relevant portions to Lu->values
// need output, turing +/- the upsampling width
// ensure that Lambda(0)>0 or Lambda(+delta)>0
void copy(Lfunc *L)
{
  bool negate_me=false; // we want f(delta)>0
  if(arb_is_negative(acb_realref(L->res[0])))
    negate_me=true;
  else
  {
    if(arb_contains_zero(acb_realref(L->res[0]))
        &&(arb_is_negative(acb_realref(L->res[1]))))
      negate_me=true;
  }

  int64_t nn=-L->u_N*L->u_stride*2; // this is where we start from
  if(negate_me) // our guess at sqrt_sign had wrong sign
    {
      acb_neg(L->sqrt_sign,L->sqrt_sign);
      for(uint64_t n=0;n<L->u_no_values;n++,nn++)
	{
	  arb_neg(L->u_values[0][n],acb_realref(L->res[nn%L->fft_NN]));
	  arb_neg(L->u_values[1][n],acb_realref(L->res[(-nn)%L->fft_NN]));
	}
    }
  else
    for(uint64_t n=0;n<L->u_no_values;n++,nn++)
    {
      arb_set(L->u_values[0][n],acb_realref(L->res[nn%L->fft_NN]));
      arb_set(L->u_values[1][n],acb_realref(L->res[(-nn)%L->fft_NN]));
    }
  for(uint64_t n=0;n<L->fft_NN;n++)
    acb_one(L->res[n]); // reclaim most of the memory
}


int64_t calc_m(uint64_t i, double two_pi_by_B, double dc)
{
  double x=log((double) i/dc)/two_pi_by_B;
  return(round(x));
}

// sks=log(m/sqrt{N})-u_m
void comp_sks(arb_t sks, uint64_t m, int64_t ms, Lfunc *L, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3;
  if(!init)
  {
    init=true;
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
  }
  arb_mul_si(tmp1,L->two_pi_by_B,ms,prec);
  arb_mul_ui(tmp2,L->one_over_root_N,m+1,prec);
  arb_log(tmp3,tmp2,prec);
  arb_sub(sks,tmp3,tmp1,prec);
}

// just check our G values go down far enough
Lerror_t finalise_comp(Lfunc *L)
{
  double two_pi_by_B=L->one_over_B*2*M_PI;
  L->offset=calc_m(1,two_pi_by_B,L->dc);
  if(L->offset<L->low_i)
  {
    // The grid floor (fixed at umin=-32 ln2, so degree-dependent only) does not
    // reach the conductor-dependent offset we need.  Report it as a recoverable
    // fatal error rather than exit()ing out from under the caller.
    if(verbose)
      fprintf(stderr,"G values need to go down to %" PRId64 ". We only have down to %" PRId64 ".\n",
          L->offset,L->low_i);
    return ERR_G_EXTENT;
  }
  return ERR_SUCCESS;
} /* finalise_comp */


// do the k convolutions, summing results into res.
void do_convolves(Lfunc *L)
{
  acb_t tmp;
  arb_t tmp1;
  acb_init(tmp);
  arb_init(tmp1);
  int64_t n, n2, prec=L->wprec;
  //int64_t n0=round(log(L->M0/Lfu->dc)/(Lf->one_over_B*2*M_PI));
  int64_t n0=L->low_i;
  if(verbose)
    printf("Taking G values from %" PRId64 " to %" PRId64 "\n",n0,L->hi_i);

  // let's do k=0
  //printf("Doing k=%" PRIu64 "\n",0);fflush(stdout);
  for(n = 0; n < (int64_t)L->fft_N; n++)
    acb_zero(L->G[n]);

  // just copy those G we actually need
  // i.e. from hi_i down to u_m=round(log(1/sqrt{conductor})*B/2/Pi)
  // forget that, copy them all
  for(n=n0,n2=n0-L->low_i;n<=L->hi_i;n++,n2++)
  {
    int64_t n1=n%L->fft_N;
    //printf("copying G value ");arb_printd(Lf->Gs[0][n2],20);printf("\n");
    arb_set(acb_realref(L->G[n1]),L->Gs[0][n2]);
    //arb_zero(acb_imagref(L->G[n1]));
  }
  acb_convolve(L->res,L->skm[0],L->G,L->fft_N,L->w,prec);
  if ( verbose ) {
    printf("Convolve 1 out of %" PRId64 " completed.\n",L->max_K);
  }
  uint64_t k;
  for(k=1;k<L->max_K;k++) {
    for(n=0;n<(int64_t)L->fft_N;n++)
      acb_zero(L->G[n]);
    for(n=n0,n2=n0-L->low_i;n<=L->hi_i;n++,n2++) {
      int64_t n1=n%L->fft_N;
      arb_set(acb_realref(L->G[n1]),L->Gs[k][n2]);
    }
    acb_convolve(L->kres,L->skm[k],L->G,L->fft_N,L->w,prec);
    if(verbose)
      printf("Convolve %" PRIu64 " out of %" PRId64 " completed.\n",k+1,L->max_K);

    for(n=0; n <= (int64_t)L->fft_N/2; n++)
      acb_add(L->res[n],L->res[n],L->kres[n],prec);
    // we need [N-1] to compute sign when F_hat(0)=0
    acb_add(L->res[L->fft_N-1],L->res[L->fft_N-1],L->kres[L->fft_N-1],prec);
  }
  acb_clear(tmp);
  arb_clear(tmp1);
} /* do_convolves */

// handle the coefficients from m=1 to M0-1
// that weren't convolved
void finish_convolves(Lfunc *L)
{
  if(L->M0==1) return;

  int64_t prec=L->wprec;
  acb_t tmp,tmp2;
  arb_t tmp1,sks;
  acb_init(tmp);acb_init(tmp2);
  arb_init(tmp1);arb_init(sks);

  if(verbose)
  {
    printf("Finishing off convolutions...\n");
    printf("F_hat[-1]=");acb_printd(L->res[L->fft_N-1],30);printf("\n");
  }
  double two_pi_by_B=2.0*M_PI*L->one_over_B;
  int64_t n;
  uint64_t m;
  for(m=0;m<L->M0-1;m++)
  {
    //printf("Doing m=%" PRIu64 "\n",m);
    int64_t ms=calc_m(m+1,two_pi_by_B,L->dc);
    for(n=-1;;n++)
    {
      int64_t nn=ms+n;
      if(nn>L->hi_i) // run out of G values
        break;
      acb_mul_arb(tmp2,L->ans[m],L->Gs[0][nn-L->low_i],prec);
      //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
      acb_add(L->res[n%L->fft_N],L->res[n%L->fft_N],tmp2,prec);
    }
    comp_sks(sks,m,ms,L,prec);
    for(n=-1;;n++)
    {
      int64_t nn=ms+n;
      if(nn>L->hi_i) // run out of G values
        break;
      acb_mul_arb(tmp,L->ans[m],sks,prec);
      acb_mul_arb(tmp2,tmp,L->Gs[1][nn-L->low_i],prec);
      //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
      acb_add(L->res[n%L->fft_N],L->res[n%L->fft_N],tmp2,prec);
    }
    arb_set(tmp1,sks);
    uint64_t k;
    for(k=2;k<L->max_K;k++)
    {
      arb_mul(tmp1,tmp1,sks,prec);
      acb_mul_arb(tmp,L->ans[m],tmp1,prec); // an/sqrt(n)(log(m/sqrt(N))-um)^k
      for(n=-1;;n++)
      {
        int64_t nn=ms+n;
        if(nn>L->hi_i) // run out of G values
          break;
        acb_mul_arb(tmp2,tmp,L->Gs[k][nn-L->low_i],prec);
        //if(n==-1) {printf("adding ");acb_printd(tmp2,20);printf("\n");}
        acb_add(L->res[n%L->fft_N],L->res[n%L->fft_N],tmp2,prec);
      }
    }
  }

  if(verbose)
    printf("Convolutions finished.\n");

  acb_clear(tmp);acb_clear(tmp2);
  arb_clear(tmp1);arb_clear(sks);
} /* finish_convolves */

// do the final iFFT (with implicit upsampling)
// multiply output by 2Pi/B
// add the ftwiddle error
void final_ifft(Lfunc *L)
{
  if(verbose)
  {
    printf("Going to do length %" PRIu64 " iFFT.\n",L->fft_NN);
    fflush(stdout);
  }
  int64_t prec=L->wprec;
  if(verbose)
  {
    printf("Pre-iFFT:-\n");
    for(uint64_t i=0;i<48;i++)
    {
      printf("res[%" PRIu64 "] = ",i);
      acb_printd(L->res[i],20);
      printf("\n");
    }
  }
  acb_ifft(L->res,L->fft_NN,L->ww,prec);
  if(verbose){printf("iFFT done.\n");fflush(stdout);}

  for(uint64_t n=0;n<L->fft_NN;n++)
  {
    arb_mul(acb_realref(L->res[n]),acb_realref(L->res[n]),L->two_pi_by_B,prec);
    arb_add_error(acb_realref(L->res[n]),L->ftwiddle_error);
  }
} /* final_ifft */

// Build M = L^(1/k), compute it, and materialize its data into L. Returns the
// error code from M's computation (fatal bits abort the assembly).
static Lerror_t extract_and_assemble(Lfunc *L, uint64_t k)
{
  // ---- M's parameters: degree/k, cond^(1/k), mus every k-th ----
  Lparams_t Mp;
  Mp.degree = L->degree / k;
  ulong rem; Mp.conductor = (uint64_t) n_rootrem(&rem, (ulong)L->conductor, (ulong)k);
  Mp.normalisation = L->normalisation;    // M shares L's algebraic->analytic shift
  double *mmus = (double *) malloc(sizeof(double)*Mp.degree);
  for (uint64_t i = 0; i < Mp.degree; i++)
    mmus[i] = L->mus[i*k] - L->normalisation; // M's algebraic mus (L->mus is analytic)
  Mp.mus = mmus;
  // CRITICAL: M is fed the RAW (algebraic) factors retained_f^(1/k), so M must carry
  // the same normalisation L was created with; M->use_lpoly then re-normalises them
  // exactly as L's did. (Setting normalisation=0 here computes M in the wrong
  // normalization and produces garbage / RH errors.)
  Mp.target_prec = L->target_prec; Mp.wprec = 0; Mp.gprec = 0;
  Mp.self_dual = L->self_dual; Mp.rank = DK; Mp.cache_dir = L->cache_dir;
  Mp.extract_powers = YES;                // inherit -> free recursion for M^(k*k')

  Lerror_t ec = ERR_SUCCESS;
  Lfunc_t Mt = Lfunc_init_advanced(&Mp, &ec);
  free(mmus);
  if (fatal_error(ec)) return ec;         // e.g. ERR_BAD_DEGREE when M would be degree 1
  Lfunc *M = (Lfunc *) Mt;

  // ---- supply M_p = retained_f^(1/k) for primes up to nmax(M) ----
  // EVERY supplied factor (good AND bad) was certified an exact integer k-th power by
  // power_extract_prepare, so re-extracting the exact integer root here is sound (a
  // false L = M^k would have already been rejected with ERR_POWER). We feed M the EXACT
  // fmpz_poly root, not a ball k-th root, so M sees the same point factors L did.
  uint64_t Mmax = Lfunc_nmax(Mt);
  acb_poly_t Mp_poly; acb_poly_init(Mp_poly);
  fmpz_poly_t Mp_z; fmpz_poly_init(Mp_z);
  for (uint64_t i = 0; i < L->n_retained; i++) {
    uint64_t p = L->retained_p[i];
    if (p > Mmax) continue;               // retention is in SUPPLY order, not prime
                                          // order (callers may append bad factors
                                          // last), so skip-and-continue, never break.
    poly_exact_kth_root(Mp_z, &L->retained_f[i], k, M->wprec); // certified above
    acb_poly_set_fmpz_poly(Mp_poly, Mp_z, M->wprec);
    Lfunc_use_lpoly(Mt, p, Mp_poly);
  }
  fmpz_poly_clear(Mp_z); acb_poly_clear(Mp_poly);

  ec |= Lfunc_compute(Mt);
  if (fatal_error(ec)) { Lfunc_clear(Mt); return ec; }

  // ---- materialize M^k into L ----
  L->rank = (int64_t) k * M->rank;
  // Emit each of M's isolated ordinates k times. find_zeros fills zeros[side][0..]
  // with positive ordinates; unfilled trailing slots stay exactly 0 (arb_is_zero),
  // which marks the end. (For self_dual M, side 1 is empty -> L side 1 is empty too.)
  for (uint64_t side = 0; side < 2; side++) {
    uint64_t dst = 0;
    for (uint64_t src = 0; src < MAX_ZEROS && dst < MAX_ZEROS; src++) {
      if (arb_is_zero(M->zeros[side][src])) break; // no more isolated zeros
      for (uint64_t r = 0; r < k && dst < MAX_ZEROS; r++)
        arb_set(L->zeros[side][dst++], M->zeros[side][src]);
    }
  }
  acb_pow_ui(L->sign, M->sign, k, L->wprec);          // eps^k
  acb_pow_ui(L->sqrt_sign, M->sqrt_sign, k, L->wprec); // sqrt_sign^k
  arb_pow_ui(L->L_d, M->L_d, k, L->wprec);            // leading Taylor coeff ^ k
  // Lam_d (= Lambda^(rank)(1/2), an unnormalized derivative) is intentionally NOT
  // maintained for an assembled L: the pipeline stage that consumes it is skipped,
  // L_d is set directly above, and Lfunc_Taylor returns L_d. (It is not a plain k-th
  // power of M's: the correct value would be (k*r)!/(r!)^k * Lam_d(M)^k.)

  // ---- attach the factor (owned by L; cleared in clear.c) ----
  L->factors = (Lfunc_t *) malloc(sizeof(Lfunc_t));
  L->factor_mults = (uint64_t *) malloc(sizeof(uint64_t));
  L->factors[0] = Mt; L->factor_mults[0] = k; L->n_factors = 1;
  return ec;                                          // may carry warnings
}

// this is called by the user to compute all the bits of the Lfunc we expect them to want
// including Lambda(t) for t =0,1/A,2/A,....
// the zeros up to height 64/degree
// the (apparent) rank
// sign and sqrt_sign
// Lambda^(rank)(1/2)
// L^(rank)(1/2)/rank!
Lerror_t Lfunc_compute(Lfunc_t Lf)
{
  static bool init=false;
  static arb_t tmp1,sks;
  static acb_t ctmp;
  if(!init)
  {
    init=true;
    arb_init(tmp1);
    arb_init(sks);
    acb_init(ctmp);
  }

  Lfunc *L=(Lfunc *) Lf;

  if(verbose) {
    for(int i = 0; i < 100; i++) {
      printf("a[%d] = ", i + 1);
      acb_printd(L->ans[i], 6);
      printf("\n");
    }
  }

  int64_t prec=L->wprec;

  // Reject perfect powers / repeated-factor L-functions up front: their doubled
  // zeros break zero-finding (stat_point -> ERR_DBL_ZERO -> RH count mismatch).
  // The caller can opt in to extraction via Lparams.extract_powers.
  {
    Lerror_t guard = power_guard(L);
    if (fatal_error(guard)) {
      if (L->extract_powers != YES)
        return guard;                    // opted out: reject as before
      uint64_t k = 0;
      Lerror_t prep = power_extract_prepare(L, &k);
      if (fatal_error(prep))
        return prep;                     // confirmed-not-a-clean-power -> ERR_POWER
      return extract_and_assemble(L, k); // build M, compute, assemble; done
    }
  }

  #ifdef BUTHE
  buthe_Wf_error(L); // add the error for the missing tail
  if(verbose){printf("Buthe Wf = ");arb_printd(L->buthe_Wf,20);printf("\n");fflush(stdout);}
  #endif

  // when we get here, the normalised L->M dirichlet coefficients are in L->ans[0]..[M-1]
  // use the first M0 of them
  arb_zero(L->sum_ans);
  for(uint64_t m=0;m<L->M0-1;m++) // we need to divide by sqrt(n) NOT a normalisation
  {
    arb_sqrt_ui(tmp1,m+1,prec); // n^(1/2)
    acb_div_arb(L->ans[m],L->ans[m],tmp1,prec);
    acb_abs(tmp1,L->ans[m],prec);
    arb_add(L->sum_ans,L->sum_ans,tmp1,prec);
  }
  if(verbose){printf("sum_{n <= %"  PRIu64 " |an/sqrt(n)|=",L->M0 -1);arb_printd(L->sum_ans,10);printf("\n");fflush(stdout);}
  for(uint64_t k=0;k<L->max_K;k++)
    for(uint64_t n=0;n<L->fft_N;n++)
      acb_zero(L->skm[k][n]);

  double two_pi_by_B=2.0*M_PI*L->one_over_B;

  for(uint64_t m=L->M0-1;m<L->M;m++)
  {
    int64_t ms=calc_m(m+1,two_pi_by_B,L->dc);
    int64_t b=(-ms)%L->fft_N;
    //arb_set_ui(n,m+1);
    arb_sqrt_ui(tmp1,m+1,prec);
    acb_div_arb(L->ans[m],L->ans[m],tmp1,prec);
    acb_abs(tmp1,L->ans[m],prec);
    arb_add(L->sum_ans,L->sum_ans,tmp1,prec);
    //printf("sum |an| now ");arb_printd(L->sum_ans,10);printf("\n");
    acb_add(L->skm[0][b],L->skm[0][b],L->ans[m],prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^0
    comp_sks(sks,m,ms,L,prec);
    //printf("m=%" PRIu64 " sks=",m);arb_printd(sks,20);printf("\n");
    acb_mul_arb(ctmp,L->ans[m],sks,prec);
    acb_add(L->skm[1][b],L->skm[1][b],ctmp,prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^1
    arb_set(tmp1,sks);
    for(uint64_t k=2;k<L->max_K;k++)
    {
      arb_mul(tmp1,tmp1,sks,prec);
      acb_mul_arb(ctmp,L->ans[m],tmp1,prec);
      acb_add(L->skm[k][b],L->skm[k][b],ctmp,prec); // a_m/sqrt(m)(log(m/sqrt(N))-u_m)^k
    }
    //if( (m + 1) % 100 == 0){printf("sum_{n <= %ld |an/sqrt(n)|=",m+1);arb_printd(L->sum_ans,10);printf("\n");fflush(stdout);}
  }
  if(verbose){printf("sum_{n <= %"  PRIu64 " |an/sqrt(n)|=",L->M);arb_printd(L->sum_ans,10);printf("\n");fflush(stdout);}
  Lerror_t ecode=finalise_comp(L);
  if(fatal_error(ecode))
    return ecode;
  do_convolves(L);
  finish_convolves(L);
  ecode|=do_pre_iFFT_errors(L);
  if(fatal_error(ecode))
    return ecode;
  final_ifft(L);

  if(verbose)
    for(uint64_t i=0;i<=L->fft_NN/OUTPUT_RATIO;i+=L->fft_NN/128)
    {
      arb_set_d(acb_realref(ctmp),0.5);
      arb_mul_ui(acb_imagref(ctmp),L->one_over_A,i,L->wprec);
      abs_gamma(tmp1,ctmp,L,L->wprec);
      arb_div(sks,acb_realref(L->res[i]),tmp1,L->wprec);
      printf("Lambda(1/2+i%f)=",(double)i/L->A);arb_printd(sks,20);printf("\n");
    }

  copy(L);

#ifdef COMPUTE_RANK
  ecode|=do_rank(L);
  if(fatal_error(ecode))
    return ecode;
#endif

#ifdef COMPUTE_ZEROS
  arb_set_d(acb_realref(ctmp),0.5);
  arb_zero(acb_imagref(ctmp));
  abs_gamma(sks,ctmp,L,L->wprec);
  arb_div(L->L_d,L->Lam_d,sks,L->wprec);
  for(int i = 2; i <= L->rank; i++)
    arb_div_ui(L->L_d,L->L_d,i,L->wprec);

  ecode|=find_zeros(L,0);
  if(fatal_error(ecode))
    return ecode;
  if(L->self_dual!=YES) // don't know or definately not
  {
    ecode|=find_zeros(L,1);
    if(fatal_error(ecode))
      return ecode;
  }
#ifdef BUTHE
  ecode|=buthe_check_RH(L);
#endif

#ifdef TURING
  ecode|=turing_check_RH(L,prec);
#endif
  
#endif

  return ecode;

}

#ifdef __cplusplus
}
#endif
