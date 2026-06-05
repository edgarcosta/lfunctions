#include <flint/acb_poly.h>
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif


  // We use Buthe's method to verify our list of zeros
  // we set h=8
  // and use the zeros up to height 96/r to check the list to height buthe_b=64/r
  // This is probably huge overkill
  //
  // Caveat: Buthe's general inequality has an extra pole contribution. This
  // implementation assumes the library's current scope of entire L-functions;
  // do not use it for objects with zeta/principal poles unless that term is
  // added to S.


  // Buthe smoothing parameter h.  Was emitted by gp/buthe_ints.gp into
  // gp/buthe_ints.out alongside the (now removed) archimedean-integral table;
  // the integrals are computed on the fly in buthe_Winf via buthe_winf_integral,
  // so we keep only this scalar here.
#define BUTHE_H (8)

  // Descending grid of Buthe smoothing parameters h. Lower h lowers the
  // completeness statistic S for high-degree/high-conductor L (sweep: deg 6-9),
  // at the cost of a larger W_f tail error (so it self-limits for small-M, low-
  // degree objects). grid[0]=8 is the historical default (deg 2-5 certify there);
  // h<2*pi*b/5 holds for all grid h at every degree 2-9 (checked in init_buthe).
  static const double buthe_h_grid[BUTHE_NH] = {8.0, 6.0, 5.0, 4.0, 3.0};

  // setup Buthe zero check stuff
  Lerror_t init_buthe(Lfunc *L, int64_t prec) {
    for(int i=0;i<BUTHE_NH;i++)
      arb_init(L->buthe_Wf[i]);
    arb_init(L->buthe_Winf);
    arb_init(L->buthe_Ws);
    arb_init(L->buthe_b); // will confirm RH in [0,b]
    arb_init(L->buthe_C);
    arb_init(L->buthe_h);


    arb_set_ui(L->buthe_C,L->degree);
    arb_div_ui(L->buthe_b,L->B,OUTPUT_RATIO,prec);
    //printf("Buthe b set to ");arb_printd(L->buthe_b,20);printf("\n");
    //
    // we set buthe_h as big as possible
    // we need b-a>5h/Pi were b-a =512/16/10
    // 
    arb_set_ui(L->buthe_h,BUTHE_H);
    //printf("Buthe h set to ");arb_printd(L->buthe_h,20);printf("\n");

    // check h<pi*2*b/5 will be OK so long as degree <=10
    arb_t tmp;
    arb_init(tmp);
    arb_mul(tmp,L->pi,L->buthe_b,prec);
    arb_div_ui(tmp,tmp,5,prec);
    arb_mul_2exp_si(tmp,tmp,1);
    if(verbose) {
      printf("Max allowed h is ");
      arb_printd(tmp,20);
      printf("\n");
      fflush(stdout);
    }
    arb_sub(tmp,tmp,L->buthe_h,prec);
    if(!arb_is_positive(tmp))
    {
      arb_clear(tmp);
      return ERR_BUTHE_PARAMS;
    }

    arb_clear(tmp);
    return ERR_SUCCESS;
  }

  // the terms making up wf are computed using the coefficients bm
  // whilst processing the Euler Factors. see use_inv_lpoly in coeff.c
  //
  // compute a term for wf
  // bm is Re([T^m] T P'(T)/P(T)), where the local factor is 1/P(T).
  // Thus bm is the negative of the usual coefficient A_m in
  // -L'/L(s)=sum_{p,m} log(p) A_m p^{-ms}. Buthe's W_f has a leading
  // minus sign, so the contribution for the symmetric interval [-b,b] is
  // +2 bm sin(b m log(p))/(pi m sqrt(p^m) cosh(h m log(p)/2)).
  void wf1(Lfunc *L, uint64_t m, uint64_t pm, arb_t bm, int64_t prec)
  {
    if(arb_is_zero(bm)) // nothing to do
      return;

    arb_t s,tmp,tmp1,tmp2,logp,logpm,num,hcur;
    arb_init(s);
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(logp);
    arb_init(logpm);
    arb_init(num);
    arb_init(hcur);


    arb_log_ui(logpm,pm,prec);

    // h-independent numerator num = 2 bm sin(b m log p) / (sqrt(p^m) pi m).
    // Each grid h then contributes num / cosh(h m log p / 2) to buthe_Wf[i];
    // grid[0]=8 reproduces the historical single-h term.
    arb_sqrt_ui(tmp,pm,prec);
    arb_mul(tmp1,tmp,L->pi,prec);
    arb_mul_ui(tmp1,tmp1,m,prec); // sqrt(p^m) pi m
    arb_mul(tmp,L->buthe_b,logpm,prec);
    arb_sin(s,tmp,prec);
    arb_mul(tmp,s,bm,prec); // b(p^m) sin(.)
    arb_div(num,tmp,tmp1,prec);
    arb_mul_2exp_si(num,num,1); // 2 bm sin(.) / (sqrt(p^m) pi m)

    for(int i=0;i<BUTHE_NH;i++)
    {
      arb_set_d(hcur,buthe_h_grid[i]); // grid h are exact integers, so set_d is exact
      arb_mul(tmp1,logpm,hcur,prec);
      arb_mul_2exp_si(tmp1,tmp1,-1);
      arb_cosh(tmp2,tmp1,prec); // cosh(h m/2 log p)
      arb_div(tmp1,num,tmp2,prec);
      arb_add(L->buthe_Wf[i],L->buthe_Wf[i],tmp1,prec);
    }

    arb_clear(s);
    arb_clear(tmp);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(logp);
    arb_clear(logpm);
    arb_clear(num);
    arb_clear(hcur);

    return;
  }


  // fp1 is f^-1, fp is f
  void wf(Lfunc *L, uint64_t p, acb_poly_t fp1, acb_poly_t fp, int64_t prec)
  {
    static bool init=false;
    static acb_t tmp1,acm;
    if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(acm);
    }

    uint64_t pm=p,m=1;
    //if(verbose) if(p==2) {printf("p=%lu\n",p);acb_poly_printd(fp,20);printf("\n--------\n");acb_poly_printd(fp1,20);printf("\n--------\n");}
    while(pm<=L->buthe_M)
    {
      acb_zero(acm);
      // this is lazy. if m-1 exceeds length of fp1, adjust i,j accordingly
      for(int64_t i=1,j=m-1;(i<acb_poly_length(fp))&&(j>=0);i++,j--)
      {
        if(j>=acb_poly_length(fp1)) continue; // and avoid this
        acb_mul(tmp1,acb_poly_get_coeff_ptr(fp1,j),acb_poly_get_coeff_ptr(fp,i),prec);
        acb_mul_si(tmp1,tmp1,i,prec);
        acb_add(acm,acm,tmp1,prec);
      }
      wf1(L,m,pm,acb_realref(acm),prec);
      pm*=p;
      m++;
    }
  }

  // We assume that L=exp(sum c(p^m) p^{-ms})
  // and that |c(p^m)|<=rp^(m/2)
  // sigma_0=1, sigma_1=3/2, C=r
  // see Lemma 3.4 of Buthe
  void buthe_Wf_error(Lfunc *L)
  {
    int64_t prec=L->wprec;
    uint64_t r=L->degree;
    arb_t logM,tmp,tmp1,tmp2,hcur;
    arb_init(logM);
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(hcur);
    if(verbose)
      printf("In buthe_Wf_error with M=%lu\n",L->buthe_M);

    arb_log_ui(logM,L->buthe_M,prec); // h-independent: computed once

    // tail bound 8 r M^((2-h)/2) / ((h-2) pi) per grid h, added to buthe_Wf[i].
    for(int i=0;i<BUTHE_NH;i++)
    {
      arb_set_d(hcur,buthe_h_grid[i]); // grid h exact integers, so set_d is exact
      if(verbose) {printf("   and h = ");arb_printd(hcur,20);printf("\n");}
      arb_sub_ui(tmp1,hcur,2,prec); // h-2 = 6
      if(verbose) {printf("h-2 = ");arb_printd(tmp1,10);printf("\n");}
      arb_mul_2exp_si(tmp1,tmp1,-1); // (h-2)/2 = 3
      arb_mul(tmp2,logM,tmp1,prec); // (h-2)/2 log M
      arb_neg(tmp2,tmp2); // (2-h)/2 log M
      arb_exp(tmp,tmp2,prec); // M^((2-h)/2)
      if(verbose)
        {
          printf("M^((2-h)/2)=");
          arb_printd(tmp,20);
          printf("\n");
        }
      arb_div(tmp2,tmp,tmp1,prec); // 2M^()/(h-2)
      arb_mul_2exp_si(tmp2,tmp2,2); // 8 M^()/(h-2)
      arb_mul_ui(tmp,tmp2,r,prec); // 8 r M^()/(h-2)
      arb_div(tmp1,tmp,L->pi,prec); // /pi
      if(verbose){printf("error in Buthe Wf <= ");arb_printd(tmp1,20);printf("\n");}
      arb_add_error(L->buthe_Wf[i],tmp1);
    }
    arb_clear(logM);
    arb_clear(tmp);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(hcur);

  }

  void buthe_fhat(arb_t res, arb_t z, Lfunc *L, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,tmp3;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
    arb_div(tmp,L->pi,L->buthe_h,prec); // Pi/h
    arb_add(tmp1,z,L->buthe_b,prec); // z+b
    arb_mul(tmp2,tmp,tmp1,prec); // Pi/h z
    arb_exp(tmp1,tmp2,prec); // exp(Pi/h (z+b))
    arb_atan(tmp2,tmp1,prec); // atan(exp(Pi/h (z+b)))
    arb_sub(tmp1,z,L->buthe_b,prec); // z-b
    arb_mul(tmp3,tmp1,tmp,prec); // Pi/h(z-b)
    arb_exp(tmp,tmp3,prec); // exp(Pi/h(z-b))
    arb_atan(tmp3,tmp,prec); // atan(exp(Pi/h(z-b)))
    arb_sub(tmp,tmp2,tmp3,prec);
    arb_div(res,tmp,L->pi,prec);
    arb_mul_2exp_si(res,res,1); // 2/Pi*(atan - atan)
  }


  // sum over zeros for a dual l-function
  void buthe_Ws_dual(arb_t res, Lfunc *Lf, arb_t *zeros, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
    }
    for(uint64_t z=0;;z++)
    {
      if((z==MAX_ZEROS)||(arb_is_zero(zeros[z])))
        break;
      arb_set(tmp1,zeros[z]);
      buthe_fhat(tmp,tmp1,Lf,prec);
      if(verbose) {
        printf("Zero at ");arb_printd(tmp1,20);printf(" contributed ");arb_printd(tmp,20);printf("\n");
      }
      arb_add(res,res,tmp,prec);
    }
    arb_mul_2exp_si(res,res,1);
    if(Lf->rank>0)
    {
      arb_zero(tmp1);// we assume the rank is correct so zeros are at exactly 0
      buthe_fhat(tmp,tmp1,Lf,prec);
      arb_mul_ui(tmp1,tmp,Lf->rank,prec);
      arb_add(res,res,tmp1,prec);
    }
  }

  // sum over zeros for a non-dual l-function (will work with dual as well)
  void buthe_Ws_non_dual(arb_t res, Lfunc *Lf, arb_t *zeros, uint64_t side, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1;
    if(!init) {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
    }
    for(uint64_t z=0;;z++)
    {
      if((z==MAX_ZEROS)||(arb_is_zero(zeros[z])))
        break;
      arb_set(tmp1,zeros[z]);
      buthe_fhat(tmp,tmp1,Lf,prec);
      if(verbose) {
        printf("Zero at ");arb_printd(tmp1,20);printf(" contributed ");arb_printd(tmp,20);printf("\n");
      }
      arb_add(res,res,tmp,prec);
    }
    if((side==0)&&(Lf->rank>0)) // there are some central zeros so include them
    {
      arb_zero(tmp1); // we assume the rank is correct so zeros are at exactly 0
      buthe_fhat(tmp,tmp1,Lf,prec);
      arb_mul_ui(tmp1,tmp,Lf->rank,prec);
      arb_add(res,res,tmp1,prec);
    }
  }

  // add digamma(1/4+mu/2)
  void buthe_lgam1(arb_t res, double mu, int64_t prec)
  {
    static bool init=false;
    static arb_t s,tmp;
    if(!init)
    {
      arb_init(s);
      arb_init(tmp);
      init=true;
    }
    arb_set_d(s,1.0/4.0+(double)mu/2.0); // all normalised to (1-s)
    arb_digamma(tmp,s,prec);
    //if(verbose) {printf("lgam1(%f) returning ",mu);arb_printd(tmp,20);printf("\n");}
    arb_add(res,res,tmp,prec);
  }

  // add log N - r log pi
  void buthe_lgam2(arb_t res, uint64_t r, uint64_t N, arb_t logpi, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,tmp3;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
    arb_log_ui(tmp1,N,prec); // log N
    arb_mul_ui(tmp2,logpi,r,prec); // log Pi^r
    arb_sub(tmp3,tmp1,tmp2,prec); // log(N/Pi^r)
    if(verbose){printf("lgam2 adding ");arb_printd(tmp3,20);printf("\n");}
    arb_add(res,res,tmp3,prec);
  }

  void buthe_Winf(arb_t res, Lfunc *L, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,logpi;
    if(!init) {
      init=true;
      arb_init(tmp);
      arb_init(logpi);
    }
    arb_log(logpi,L->pi,prec);
    arb_zero(res);
    for(uint64_t k=0;k<L->degree;k++)
      buthe_lgam1(res,L->mus[k],prec);
    buthe_lgam2(res,L->degree,L->conductor,logpi,prec);
    arb_div(tmp,L->buthe_b,L->pi,prec);
    arb_mul(res,res,tmp,prec);
    if(verbose) {printf("Winf before integral = ");arb_printd(res,20);printf("\n");}
    // archimedean integrals: one per Gamma_R factor, computed rigorously on the
    // fly for the live (b,h) (replaces the old gp/buthe_ints.out static table).
    for(uint64_t k=0;k<L->degree;k++) {
      buthe_winf_integral(tmp,L->buthe_b,L->buthe_h,L->mus[k],prec);
      if(verbose){printf("Integral contributing ");arb_printd(tmp,20);printf("\n");}
      arb_add(res,res,tmp,prec);
    }
    if(verbose) {printf("Winf = ");arb_printd(res,20);printf("\n");}
  }

  // S = Wf[i] + Winf - Ws at grid point i: sets L->buthe_h = buthe_h_grid[i] and
  // recomputes the h-dependent Ws and Winf. Exposed so the threshold test can
  // probe S per grid h.
  void buthe_S_at(arb_t S, Lfunc *L, int i, int64_t prec) {
    arb_set_d(L->buthe_h, buthe_h_grid[i]);
    arb_zero(L->buthe_Ws);
    if (L->self_dual == YES) buthe_Ws_dual(L->buthe_Ws, L, L->zeros[0], prec);
    else { buthe_Ws_non_dual(L->buthe_Ws, L, L->zeros[0], 0, prec);
           buthe_Ws_non_dual(L->buthe_Ws, L, L->zeros[1], 1, prec); }
    buthe_Winf(L->buthe_Winf, L, prec);
    arb_add(S, L->buthe_Wf[i], L->buthe_Winf, prec);
    arb_sub(S, S, L->buthe_Ws, prec);
  }

  Lerror_t buthe_check_RH(Lfunc *L)
  {
    static bool init=false;
    static arb_t two_zeros,one_zero;
    if(!init) {
      init=true;
      arb_init(two_zeros);
      arb_set_ui(two_zeros,98);
      arb_div_ui(two_zeros,two_zeros,100,100); // 0.98 pair-threshold (self-dual)
      arb_init(one_zero);
      arb_set_ui(one_zero,49);
      arb_div_ui(one_zero,one_zero,100,100);   // 0.49 single-zero threshold (self_dual != YES)
    }
    // A missed conjugate pair contributes > 0.98 (self-dual == YES, paired zeros);
    // a single missed zero contributes > 0.49. When self_dual is NO or DK,
    // zeros are not guaranteed paired, so the tighter 0.49 bar is required or a
    // single miss would false-confirm.
    arb_srcptr threshold = (L->self_dual == YES) ? two_zeros : one_zero;
    int64_t prec=L->wprec;
    if(verbose)
    {printf("Going to use Weil-Barner to confirm list of zeros.\n");fflush(stdout);}

    // Adaptive sweep over the descending smoothing grid buthe_h_grid[0..BUTHE_NH-1].
    // S(h) = Wf(h) + Winf(h) - Ws(h); lowering h sharpens the test function and
    // lowers S for a COMPLETE zero list, while a genuine miss keeps S >= threshold
    // at every grid h (a missed zero's contribution rises toward 0.5 as h drops).
    // We certify at the HIGHEST grid h that achieves upper(S) < threshold (that h
    // has the smallest W_f tail error). The hard over-count test S < 0 is judged
    // ONLY at i==0 (h=8, smallest W_f error) so the wider low-h error cannot
    // manufacture a false over-count.
    arb_t S, Sdiff;
    arb_init(S); arb_init(Sdiff);
    Lerror_t verdict = ERR_RH_ERROR; // default: no grid h certified
    for (int i = 0; i < BUTHE_NH; i++) {
      buthe_S_at(S, L, i, L->wprec);
      if(verbose)
      {
        printf("h = ");arb_printd(L->buthe_h,20);
        printf("\nWs = ");arb_printd(L->buthe_Ws,20);
        printf("\nWinf = ");arb_printd(L->buthe_Winf,20);
        printf("\nWf = ");arb_printd(L->buthe_Wf[i],20);
        printf("\nS = ");arb_printd(S,20);printf("\n");
      }
      if (i == 0 && arb_is_negative(S)) { // hard over-count at the safest h
        if(verbose) printf("Error in Weil-Barner check. Winf+Wf-Ws* must allow >=0. RH not confirmed.\n");
        verdict = ERR_BUT_ERROR; break;
      }
      arb_sub(Sdiff, S, threshold, prec);
      if (arb_is_negative(Sdiff)) { // upper(S) < threshold => certified at this (highest) h
        verdict = ERR_SUCCESS; break;
      }
    }
    if(verbose && (verdict == ERR_RH_ERROR)) printf("Looks like we've missed some zero(s).\n");
    arb_clear(S); arb_clear(Sdiff);
    return verdict;
  }

#ifdef __cplusplus
}
#endif
