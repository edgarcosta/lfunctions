#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif
#ifdef TURING
// int im loggamma(sigma+It) t=a..b
bool imint(arb_t res, arb_t sigma, arb_t a, arb_t b, uint64_t prec)
{
  if(verbose)
    {printf("Integrating with sig=");arb_printd(sigma,10);printf(" a=");arb_printd(a,10);printf(" and b=");arb_printd(b,10);printf("\n");}
  static arb_t a2,b2,s2,b2s2,a2s2,tmp1,tmp2,tmp3,shalf;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(a2);
      arb_init(b2);
      arb_init(s2);
      arb_init(b2s2);
      arb_init(a2s2);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(shalf);
    }
  arb_set_d(tmp1,-0.5);
  arb_add(shalf,sigma,tmp1,prec); // sigma-1/2

  arb_mul(s2,sigma,sigma,prec);
  arb_mul(a2,a,a,prec);
  arb_mul(b2,b,b,prec);
  arb_add(tmp1,a2,s2,prec);
  arb_log(a2s2,tmp1,prec);
  arb_add(tmp1,b2,s2,prec);
  arb_log(b2s2,tmp1,prec);
  arb_sub(tmp1,s2,sigma,prec);
  arb_sub(tmp2,tmp1,b2,prec);
  arb_mul(tmp3,tmp2,b2s2,prec); // log(b^2+sig^2)(sig^2-b^2-sig)
  arb_sub(tmp2,tmp1,a2,prec);
  arb_mul(tmp1,tmp2,a2s2,prec); // log(a^2+sig^2)(sig^2-a^2-sig)
  arb_sub(tmp2,tmp3,tmp1,prec);
  arb_sub(tmp1,b2,a2,prec);
  arb_mul_ui(tmp3,tmp1,3,prec);
  arb_add(res,tmp2,tmp3,prec);
  arb_mul_2exp_si(res,res,-2);
  if(verbose)
    {printf("value 1 = ");arb_printd(res,10);printf("\n");}

  arb_div(tmp1,b,sigma,prec);
  arb_atan(tmp2,tmp1,prec);
  arb_mul(tmp1,tmp2,b,prec); // b atan(b/sig)
  arb_div(tmp2,a,sigma,prec);
  arb_atan(tmp3,tmp2,prec);
  arb_mul(tmp2,tmp3,a,prec); // a atan(a/sig)
  arb_sub(tmp3,tmp2,tmp1,prec); // a atan(a/sigma)-b atan(b/sigma)
  arb_mul(tmp2,tmp3,shalf,prec);
  if(verbose)
    {printf("value 2 = ");arb_printd(tmp2,10);printf("\n");}

  arb_add(res,res,tmp2,prec);
  arb_neg(res,res);
  arb_div(tmp1,b,a,prec);
  arb_log(tmp2,tmp1,prec);
  arb_mul_2exp_si(tmp2,tmp2,-3);
  arb_add_error(res,tmp2);
  if(verbose)
    {printf("im_int returning ");arb_printd(res,10);printf("\n");}
  return(true);
}

// Q(s) is analytic conductor defined pg 387 col 2
  void logQ(arb_t res, acb_t s, Lfunc *L, int64_t prec)
{
  //printf("In logQ with s=");acb_printd(s,10);printf("\n");
  static bool init=false;
  static arb_t two_pi,tmp1,tmp2;
  static acb_t stmp1,stmp2,stmp3;
  if(!init)
    {
      init=true;
      arb_init(two_pi);arb_init(tmp1);arb_init(tmp2);
      acb_init(stmp1);acb_init(stmp2);acb_init(stmp3);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }
  acb_set_ui(stmp1,L->conductor);
  arb_set(acb_imagref(stmp2),acb_imagref(s));
  for(uint64_t j=0;j<L->degree;j++)
    {
      arb_set_d(tmp2,L->mus[j]);
      arb_add(acb_realref(stmp2),acb_realref(s),tmp2,prec);
      acb_mul(stmp3,stmp1,stmp2,prec);
      acb_div_arb(stmp1,stmp3,two_pi,prec);
    }
  if(verbose)
    {printf("Analytic conductor = ");acb_printd(stmp1,10);printf("\n");}
  acb_abs(tmp1,stmp1,prec);
  arb_log(res,tmp1,prec);
  if(verbose)
    {printf("LogQ returning ");arb_printd(res,10);printf("\n");}
}

uint64_t set_X(uint64_t r,double *mus,double one_over_B)
{
  double max_mu=mus[0];
  for(uint64_t j=1;j<r;j++)
    if(mus[j]>max_mu)
      max_mu=mus[j];
  printf("max mu_j = %10.8e\n",max_mu);
  max_mu+=2.5;
  double max_t=1.0/(OUTPUT_RATIO*one_over_B);
  printf("t = %10.8e\n",max_t);
  double X2=max_t*max_t-max_mu*max_mu;
  if(X2<=25.0)
    {
      printf("Error setting X (eq 4-10). A mu_j was too large. Exiting.\n");
      exit(0);
    }
  return(sqrt(X2));
}


// return upb for 2\int\limits_{t_0}^{t0+h} S(t) \dif t
// see Th 4.6
bool St_int(arb_t res, arb_t h, arb_t t0, Lfunc* L, int64_t prec)
{
  static bool init=false;
  static acb_t s;
  static arb_t Q1,Q2,l2_half,pi,rc_theta_etc,tmp,tmp1;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      acb_init(s);
      arb_init(pi);
      arb_set_d(acb_realref(s),1.5);
      arb_init(Q1);
      arb_init(Q2);
      arb_init(l2_half);
      arb_log_ui(Q1,2,prec);
      arb_set_d(Q2,-0.5);
      arb_add(l2_half,Q1,Q2,prec);
      arb_init(rc_theta_etc);
      arb_set_d(Q2,5.65056); // c_theta < 5.65055 in ARB Th 4.6
      arb_sqrt_ui(Q1,2,prec);
      arb_mul_ui(pi,Q1,L->X-5,prec);
      arb_inv(Q1,pi,prec); // 1.0/(sqrt(2)(X-5)
      arb_add(pi,Q1,Q2,prec); 
      arb_mul_ui(rc_theta_etc,pi,L->degree,prec); // c_\theta r+r/(\sqrt{2}(X-5))
      if(verbose)
	{printf("c theta bit = ");arb_printd(rc_theta_etc,10);printf("\n");}
      arb_const_pi(pi,prec);
      L->X=set_X(L->degree,L->mus,L->one_over_B);
    }
  arb_set(acb_imagref(s),t0); // 3/2+it0
  logQ(Q2,s,L,prec);
  arb_mul(tmp,Q2,l2_half,prec);
  arb_add(acb_imagref(s),acb_imagref(s),h,prec); // 3/2+it1
  logQ(Q1,s,L,prec);
  arb_mul_2exp_si(Q1,Q1,-2);
  arb_add(tmp1,tmp,Q1,prec);
  arb_add(tmp,tmp1,rc_theta_etc,prec);
  arb_div(res,tmp,pi,prec);
  arb_mul_2exp_si(res,res,1);
  return true;
}

bool turing_int(arb_t res, arb_t h, Lfunc *L, uint64_t side, int64_t prec)
{
  static bool init=false;
  static arb_t tmp,ptr,err;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(ptr);
      arb_init(err);
      arb_zero(err);
      arb_mul_2exp_si(tmp,L->one_over_A,-1);
      arb_add_error(err,tmp);
    }

  arb_zero(res);
  /*

  uint64_t n=1,nn=Ls->N*Ls->stride+Lc->NN/OUTPUT_RATIO;
  uint8_t this_sign,last_sign=sign(Ls->values[side][nn]);
  nn++;
  while(true)
    {
      this_sign=sign(Ls->values[side][nn]);
      if(this_sign!=last_sign)
	{
	  //printf("Turing Zero found between %lu and %lu.\n",n-1,n);
	  arb_set_d(tmp,(double)n-0.5);
	  arb_mul(ptr,Lc->one_over_A,tmp,prec);
	  arb_add(tmp,ptr,err,prec);
	  arb_sub(ptr,h,tmp,prec);
	  arb_add(res,res,ptr,prec);
	  last_sign=this_sign;
	}
      n++;nn++;
      if(n==Lc->NN/TURING_RATIO)
	return(true);
    }
  */
  return true;
}


bool turing_count(arb_t res, Lfunc *L, uint64_t prec)
{
  static bool init=false;
  static arb_t tint,tmp,tmp1,h,t0,sint,pi,rlogpi,t0hbit;
  if(!init)
    {
      init=true;
      arb_init(tint);
      arb_init(sint);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(h);
      arb_init(t0);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_div_ui(h,L->B,TURING_RATIO,prec);
      printf("Turing h set to ");arb_printd(h,10);printf("\n");
      arb_div_ui(t0,L->B,OUTPUT_RATIO,prec);
      printf("Turing t0 set to ");arb_printd(t0,10);printf("\n");
      arb_init(rlogpi);
      arb_log(tmp,pi,prec);
      arb_mul_ui(rlogpi,tmp,L->degree,prec);
      arb_init(t0hbit);
      arb_mul(tmp,t0,h,prec);
      arb_mul_2exp_si(tmp,tmp,1);
      arb_mul(tmp1,h,h,prec);
      arb_add(sint,tmp1,tmp,prec);
      arb_div(t0hbit,sint,pi,prec); 
      arb_mul_2exp_si(t0hbit,t0hbit,-1); // (2t0h+h^2)/2Pi
    }
  if(!turing_int(tint,h,L,0,prec))
    return(0);
  if(verbose)
    {printf("Turing int [0] = ");arb_printd(tint,20);printf("\n");}
  if(L->self_dual)
    arb_mul_2exp_si(tint,tint,1);
  else
    {
      if(!turing_int(tmp,h,L,1,prec))
	return(0);
      arb_add(tint,tint,tmp,prec);
    }
  if(verbose)
    {printf("Turing int [*] = ");arb_printd(tint,20);printf("\n");}
  if(!St_int(sint,h,t0,L,prec))
    return(0);
  if(verbose)
    {printf("St_int returned ");arb_printd(sint,10);printf("\n");}
  arb_div(tmp1,L->imint,pi,prec);
  arb_add(tmp,tmp1,sint,prec);
  arb_sub(sint,tmp,tint,prec); //  

  arb_log_ui(tmp1,L->conductor,prec);
  arb_sub(tmp,tmp1,rlogpi,prec);
  arb_mul(tmp1,tmp,t0hbit,prec);

  arb_add(tmp,tmp1,sint,prec);
  arb_div(res,tmp,h,prec);

  if(verbose)
    {printf("turing_count will try to return ");arb_printd(res,10);printf("\n");}
  return true;
}

Lerror_t turing_check_RH(Lfunc *L)
{
  return ERR_RH_ERROR;
}

#ifdef __cplusplus
}
#endif

#endif
