#include "glfunc.h"
#include "glfunc_internals.h"


#undef verbose
#define verbose false
#ifdef __cplusplus
extern "C"{
#endif

  // find the first point in the left tail for upsampling
  int64_t s_left_n(acb_ptr z, arb_t A, uint64_t prec) {
    static arb_t tmp,tmp1;
    static fmpz_t fmpz_tmp;
    static bool init=false;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      fmpz_init(fmpz_tmp);
    }
    arb_mul(tmp,acb_realref(z),A,prec);
    arb_get_mid_arb(tmp1,tmp);
    arb_floor(tmp,tmp1,prec);
    if(arb_get_unique_fmpz(fmpz_tmp,tmp)==0)
      return BAD_64;
    return fmpz_get_si(fmpz_tmp); // this is going to be 1st point in left tail
  }


  // sinc x = sin(pi x)/(pi x)
  // on entry sin_x = sin(pi x)
  Lerror_t s_sinc(acb_t res, acb_t pi_x, acb_t sin_pi_x, int64_t prec)
  {
    if(acb_contains_zero(pi_x))
      return ERR_SPEC_VALUE;
    acb_div(res,sin_pi_x,pi_x,prec);
    return ERR_SUCCESS;
  }

  // compute W(k/A)
  void W_k_A(arb_ptr res, Lfunc *L, int64_t k, int64_t prec, arb_t t0, arb_t pi_by_H2, arb_t A) {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,ka;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(ka);
    }
    //printf("t0 = ");arb_printd(t0,20);printf("\n");
    arb_set_si(tmp,k);
    // ka = // k/A
    arb_div(ka,tmp,A,prec);
    arb_sub(tmp,ka,t0,prec);
    arb_mul(tmp1,tmp,tmp,prec);
    arb_mul(tmp,tmp1,pi_by_H2,prec); // -pi (k/A-t0)^2/h^2
    arb_mul_ui(tmp1,ka,L->degree,prec);
    arb_mul(tmp2,tmp1,L->pi,prec);
    arb_mul_2exp_si(tmp2,tmp2,-2); // pi r k/(4A)
    arb_add(tmp1,tmp,tmp2,prec);
    // tmp = exp((k pi r)/(4 A) - (pi (k/A - t0)^2)/h^2)
    arb_exp(tmp,tmp1,prec);
    //printf("factor = ");arb_printd(tmp,20);printf("\n");
    if(k>=0)
      arb_mul(res,L->u_values_off[0][k],tmp,prec); // *Lambda(k/A)
    else
      arb_mul(res,L->u_values_off[1][-k],tmp,prec);
    //printf("W(%" PRId64 "/A) = W(",k);arb_printd(ka,20);printf(") = ");;arb_printd(res,20);printf("\n");
  }

  // W(n/A)sinc(Pi(Az-n))
  Lerror_t s_do_point(acb_ptr res, Lfunc *L, int64_t k, acb_t z, acb_t delta, acb_t sin_delta, arb_t A, int64_t prec, arb_t pi_by_H2)
  {
    static acb_t tmp,tmp1;
    static arb_t tmp2;
    static bool init=false;
    if(!init)
    {
      init=true;
      acb_init(tmp);
      acb_init(tmp1);
      arb_init(tmp2);
    }

    Lerror_t ecode=s_sinc(tmp,delta,sin_delta,prec);
    if(fatal_error(ecode))
      return ecode;
    //printf("s_sinc returned ");acb_printd(tmp,20);printf("\n");
    W_k_A(tmp2,L,k,prec,acb_realref(z),pi_by_H2,A);
    acb_mul_arb(res,tmp,tmp2,prec);
    return ecode;
    //printf("W(%" PRId64 "/A)sinc(Pi(Az-(%" PRId64 "))) = ",k,k);acb_printd(res,20);printf("\n");
  }

  // estimate W(1/2+iz) by upsampling off Lu->u_values_off
  Lerror_t s_upsample_stride(acb_ptr res, acb_ptr z, Lfunc *L, int64_t prec, arb_t pi_by_H2, arb_t err, int64_t N, uint64_t stride)
  {
    static arb_t A; // the A for upsampling = usual A / stride
    static acb_t diff,this_diff,term,sin_diff,neg_sin_diff;
    static bool init=false;
    if(!init)
    {
      init=true;
      arb_init(A);
      acb_init(diff);
      acb_init(this_diff);
      acb_init(term);
      acb_init(sin_diff);
      acb_init(neg_sin_diff);
    }
    //printf("-Pi/h^2 = ");arb_printd(pi_by_H2,20);printf("\n");
    arb_div_ui(A,L->arb_A,stride,prec);
    //printf("Special point upsampling A set to ");arb_printd(A,20);printf("\n");
    int64_t n=s_left_n(z,L->arb_A,prec);
    if(n==BAD_64)
      return ERR_SPEC_VALUE;
    acb_mul_arb(diff,z,A,prec);
    acb_sub_ui(diff,diff,n,prec);
    acb_mul_arb(diff,diff,L->pi,prec);
    acb_sin(sin_diff,diff,prec);
    acb_neg(neg_sin_diff,sin_diff);
    //printf("Nearest n = %" PRId64 "\n",n);
    //printf("Pi(Az-k0) = ");acb_printd(diff,20);printf("\n");
    int64_t nn=n; // offset into values

    // do nearest point
    Lerror_t ecode=s_do_point(res,L,n,z,diff,sin_diff,L->arb_A,prec,pi_by_H2);
    if(fatal_error(ecode))
      return ecode;
    if(verbose){printf("Nearest point contributed ");acb_printd(res,20);printf("\n");}
    acb_set(this_diff,diff);
    // do Lu->N points to left of left
    for(uint64_t count=0; count < (uint64_t)N; count++)
    {
      acb_add_arb(this_diff,this_diff,L->pi,prec);
      nn-=stride;
      if(count&1)
        ecode |= s_do_point(term,L,nn,z,this_diff,sin_diff,L->arb_A,prec,pi_by_H2);
      else
        ecode |= s_do_point(term,L,nn,z,this_diff,neg_sin_diff,L->arb_A,prec,pi_by_H2);
      if(fatal_error(ecode))
        return ecode;
      acb_add(res,res,term,prec);
    }

    acb_set(this_diff,diff);
    nn=n;

    // do N points to right of left
    for(uint64_t count = 0; count < (uint64_t)N; count++)
    {
      acb_sub_arb(this_diff,this_diff,L->pi,prec);
      nn+=stride;
      if(count&1)
        ecode|=s_do_point(term,L,nn,z,this_diff,sin_diff,L->arb_A,prec,pi_by_H2);
      else
        ecode|=s_do_point(term,L,nn,z,this_diff,neg_sin_diff,L->arb_A,prec,pi_by_H2);
      if(fatal_error(ecode))
        return ecode;
      acb_add(res,res,term,prec);
    }
    arb_add_error(acb_realref(res),err);
    arb_add_error(acb_imagref(res),err);
    return ecode;
  }

  //1/gamma_r(s)
  void acb_rgamma_r(acb_t res, acb_t s, Lfunc *L, uint64_t prec)
  {
    acb_t s_by_2,lg,ctmp,ctmp1;
    arb_t log_pi;
    acb_init(s_by_2);
    acb_init(ctmp);
    acb_init(lg);
    acb_init(ctmp1);
    arb_init(log_pi);
    arb_log(log_pi,L->pi,prec);
    //printf("in abs_gamma_r with s = ");acb_printd(s,10);printf("\n");
    acb_mul_2exp_si(s_by_2,s,-1); // s/2

    acb_mul_arb(ctmp,s_by_2,log_pi,prec); 
    acb_exp(ctmp1,ctmp,prec);
    acb_rgamma(ctmp,s_by_2,prec);
    acb_mul(res,ctmp,ctmp1,prec);
    //printf("gamma_r returning ");acb_printd(res,20);printf("\n");
    acb_clear(s_by_2);
    acb_clear(lg);
    acb_clear(ctmp);
    acb_clear(ctmp1);
    arb_clear(log_pi);
  }


  // gamma(s) per l.pdf with epsilon and N^1/2(s-1/2)
  void spec_rgamma(acb_t res, acb_t s, Lfunc *L, int64_t prec)
  {
    acb_t tmp1,tmp2,tmp3;
    acb_init(tmp1);
    acb_init(tmp2);
    acb_init(tmp3);
    arb_t tmp;
    arb_init(tmp);

    acb_set_ui(res,1);
    for(uint64_t j=0;j<L->degree;j++) {
      arb_set_d(tmp,L->mus[j]);
      acb_add_arb(tmp1,s,tmp,prec);
      acb_rgamma_r(tmp2,tmp1,L,prec);
      acb_mul(res,res,tmp2,prec);
    }
    arb_set_d(tmp,0.5);
    acb_sub_arb(tmp1,s,tmp,prec); // s-1/2
    acb_mul_2exp_si(tmp1,tmp1,-1); // (s-1/2)/2
    arb_log_ui(tmp,L->conductor,prec);
    acb_mul_arb(tmp2,tmp1,tmp,prec);
    acb_neg(tmp2,tmp2);
    acb_exp(tmp1,tmp2,prec);
    //printf("N^(-(s-1/2)/2) = ");acb_printd(tmp1,20);printf("\n");
    acb_mul(res,res,tmp1,prec);
    acb_mul(res,res,L->epsilon,prec);
    //printf("spec_gamma returning ");acb_printd(res,10);printf("\n");
    acb_clear(tmp1);
    acb_clear(tmp2);
    acb_clear(tmp3);
    arb_clear(tmp);
  }

  // compute L(s) into res where s is given in algebraic normalisation
  Lerror_t Lfunc_special_value(acb_t res, Lfunc_t LL, double alg_res, double alg_ims)
  {
    Lerror_t ecode=ERR_SUCCESS;
    Lfunc *L=(Lfunc *) LL;
    int64_t prec=L->wprec;
    //printf("Algebraic s = %f + i%f\n",alg_res,alg_ims);

    arb_t arb_err,tmp;
    acb_t an_s; // analytic version of s
    acb_init(an_s);
    arb_init(arb_err);
    arb_init(tmp);
    arb_set_d(tmp, L->normalisation);
    arb_set_d(arb_err, alg_res);
    arb_sub(acb_realref(an_s),arb_err,tmp,prec); // an re(s) = alg re(s)- norm
    arb_set_d(acb_imagref(an_s),alg_ims);
    //printf("Analytic s = ");acb_printd(an_s,20);printf("\n");

    acb_t z;
    acb_init(z);
    arb_set(acb_realref(z),acb_imagref(an_s));
    arb_set_d(tmp,0.5);
    arb_sub(acb_imagref(z),tmp,acb_realref(an_s),prec);
    if(verbose) {printf("special value z = ");acb_printd(z,20);printf("\n");}

    double T=512.0/L->degree/OUTPUT_RATIO; // largest T for which we have data
    uint64_t stride=32; // should be computed dynamically
    double A=L->A/(double)stride;

    double h=sqrt(1.0/A)*1.001,best_h=h;
    double H=ceil(A*A*h*h/2.0),best_H=H;
    double M=H/A;
    double dimz=0.5-(alg_res - L->normalisation); // double approx to Im z
    if(verbose) printf("A = %f Im z = %f\n",A,dimz);
    ecode|=arb_upsampling_error(arb_err,M,H,h,A,L->mus,L->degree,L->conductor,T,acb_imagref(z),0,L->pi,prec); // first effort
    if(fatal_error(ecode))
    {
      arb_clear(tmp);
      acb_clear(z);
      acb_clear(an_s);
      arb_clear(arb_err);
      return ecode;
    }

    arb_t best_err;
    arb_init(best_err);
    arb_set(best_err,arb_err);
    while(true)
    {
      int64_t extra_bits=(int64_t)(M_PI*(dimz*dimz/(h*h)+fabs(dimz)*A)/M_LN2)+10;
      if(verbose)
      {
        printf("h=%f extra_bits=%" PRId64 "\n", h, extra_bits);
        printf("arb_error now ");arb_printd(arb_err,20);printf("\n");
      }
      arb_mul_2exp_si(tmp,arb_err,L->target_prec+extra_bits);
      arb_sub_ui(tmp,tmp,1,prec);
      if(arb_is_negative(tmp)) // achieved target error
        break;
      h*=1.01;
      H=ceil(A*A*h*h/2.0);
      if(H*stride>L->u_no_values_off) // run out of data
	{
	  arb_set(arb_err,best_err);
	  h=best_h;H=best_H;
	  //ecode|=ERR_SPEC_PREC; // not necessarily, wait till end
	  break;
	}

      M=H/A;
      ecode|=arb_upsampling_error(arb_err,M,H,h,A,L->mus,L->degree,L->conductor,T,acb_imagref(z),0,L->pi,prec);
      if(fatal_error(ecode))
      {
        arb_clear(best_err);
        arb_clear(tmp);
        acb_clear(z);
        acb_clear(an_s);
        arb_clear(arb_err);
        return ecode;
      }
      arb_sub(tmp,best_err,arb_err,prec);
      if(arb_is_positive(tmp)) // found a better h,H so use them
      {
        best_h=h;
        best_H=H;
        arb_set(best_err,arb_err);
      }
      if(verbose)
	{
	  printf("Upsampling error now ");arb_printd(arb_err,20);printf("\n");
	}
    }
    arb_clear(best_err);

    uint64_t iH=H;
    if(verbose) printf("H = %" PRIu64 " h = %f\n",iH,h);
    if(verbose) {printf("Upsample error set to ");arb_printd(arb_err,20);printf("\n");}

    arb_t pi_by_H2;
    arb_init(pi_by_H2);
    arb_set_d(tmp,h);
    arb_mul(tmp,tmp,tmp,prec);
    arb_div(pi_by_H2,L->pi,tmp,prec);
    arb_neg(pi_by_H2,pi_by_H2); // -Pi/h^2
    ecode|=s_upsample_stride(res, z, L, prec, pi_by_H2, arb_err, iH, stride);
    if(fatal_error(ecode))
    {
      arb_clear(pi_by_H2);
      arb_clear(tmp);
      acb_clear(z);
      acb_clear(an_s);
      arb_clear(arb_err);
      return ecode;
    }
    if(verbose){printf("W(z) = ");acb_printd(res,20);printf("\n");}
    // go from W(z)->Lam(1/2+iz)
    acb_t s,ctmp;
    acb_init(s);acb_init(ctmp);
    acb_mul_arb(s,z,L->pi,prec);
    acb_mul_ui(ctmp,s,L->degree,prec);
    acb_mul_2exp_si(ctmp,ctmp,-2); // pi r z /4
    arb_mul(tmp,acb_imagref(z),acb_imagref(z),prec);
    arb_mul(tmp,tmp,pi_by_H2,prec);
    arb_sub(acb_realref(ctmp),acb_realref(ctmp),tmp,prec);
    acb_exp(s,ctmp,prec);
    if(verbose){printf("going from W to Lambda by dividing by ");acb_printd(s,20);printf("\n");}
    acb_div(res,res,s,prec);
    // go from Lam->L by dividing out gamma(s)
    spec_rgamma(ctmp,an_s,L,prec);
    if(verbose){printf("going from Lambda to L by multiplying by ");acb_printd(ctmp,20);printf("\n");}
    acb_mul(res,res,ctmp,prec);
    acb_abs(tmp,res,prec);
    arb_get_rad_arb(tmp,tmp);
    arb_mul_2exp_si(tmp,tmp,L->target_prec);
    arb_sub_ui(tmp,tmp,1,prec);
    if(!arb_is_negative(tmp))
      ecode|=ERR_SPEC_PREC;

    acb_clear(s);
    acb_clear(ctmp);
    acb_clear(z);
    acb_clear(an_s);
    arb_clear(pi_by_H2);
    arb_clear(tmp);
    arb_clear(arb_err);

    return ecode;
  }

#ifdef __cplusplus
}
#endif
