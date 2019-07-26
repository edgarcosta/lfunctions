#include "glfunc.h"
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif

  /*
  // see ARB
  uint64_t guess_rank(Lfunc *L, uint64_t side, uint64_t prec) 
  {
  static bool init=false;
  static arb_t t;
  double ratio;

  if (!init) 
  {
  init=true;
  arb_init(t);
  }
  arb_div(t,L->u_values_off[side][2],L->u_values_off[side][1],prec);
  ratio = arf_get_d(arb_midref(t),ARF_RND_NEAR);
  if(verbose)
  {
  printf("Ratio = %10.8e\n",ratio);
  printf("eps^2 = ");arb_printd(acb_realref(L->epsilon_sqr),20);printf("\n");
  }

  if(L->self_dual==NO)
  {
  if (ratio <= 1) return 0;
  return (int)floor(log(ratio)/log(2.0)+0.5);
  } 

  if (arb_is_negative(acb_realref(L->epsilon_sqr)))
  {
  if (ratio <= 2) return 1;
  return 1+2*(int)floor(log(ratio/2)/log(4.0)+0.5);
  } 
  else 
  {
  if (ratio <= 1) return 0;
  return 2*(int)floor(log(ratio)/log(4.0)+0.5);
  }
  }
  */

  // exp(Pi(r*n/(4*A)-n^2/(A^2*H^2)))
  void exp_term(arb_t res,int64_t n,Lfunc *L,int64_t prec)
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
    arb_mul_si(tmp1,L->u_one_over_A,n,prec); // n/A
    arb_mul(tmp2,tmp1,tmp1,prec); // n^2/A^2
    arb_mul(tmp3,tmp2,L->u_pi_by_H2,prec); // -Pi n^2/(A^2H^2)
    arb_mul_ui(tmp2,tmp1,L->degree,prec); //  r n/A
    arb_mul(tmp1,tmp2,L->pi,prec);  // Pi r n/A
    arb_mul_2exp_si(tmp1,tmp1,-2); // Pi r n/(4A)
    arb_add(tmp2,tmp1,tmp3,prec);
    arb_exp(res,tmp2,prec);
  }


  // compute Lambda^(d)(1/2)
  // assumes Lambda^(k)(1/2)=0 for all k in [0,d-1]
  bool df_zero(arb_t res, uint64_t d, Lfunc *L, int64_t prec) {
    static int64_t init_prec=0;
    static bool init=false;
    static arb_t a_pi_d[MAX_L+1], d_bang_pi[MAX_L+1], tmp, tmp1, an,term;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(term);
      arb_init(an);
      for(uint64_t i=0;i<=MAX_L;i++)
      {arb_init(a_pi_d[i]);arb_init(d_bang_pi[i]);}
    }
    if(prec>init_prec) // precision has increased
    {
      init_prec=prec;
      arb_set_ui(a_pi_d[0],1);
      arb_set(a_pi_d[1],L->u_pi_A);
      for(uint64_t i=2;i<=MAX_L;i++)
        arb_mul(a_pi_d[i],a_pi_d[i-1],L->u_pi_A,prec);
      arb_set_ui(d_bang_pi[0],1);
      arb_inv(d_bang_pi[1],L->pi,prec); // d!/Pi^d
      for(uint64_t i=2;i<=MAX_L;i++)
      {
        arb_mul_ui(d_bang_pi[i],d_bang_pi[i-1],i,prec);
        arb_div(d_bang_pi[i],d_bang_pi[i],L->pi,prec);
      }
    }

    if(d>MAX_L) {
      printf("Can't compute F^(d)(0) for d>%d.\n",MAX_L);
      return false;
    }
    if(d==0)
    {
      arb_set(res,L->u_values_off[0][0]);
      return true;
    }
    // do the n=0 term. We know L(1/2)=0
    arb_zero(res);
    for(uint64_t n = 1; n <= L->u_N; n++)
    {
      arb_set_ui(an,n);
      arb_zero(term);
      bool neg_me;
      if(n&1)
        neg_me=false;
      else
        neg_me=true;
      for(int64_t D=d;D>0;D-=2)
      {
        arb_pow_ui(tmp,an,D,prec);
        arb_div(tmp1,d_bang_pi[D],tmp,prec); // D!/(Pi n)^D
        if(neg_me)
        {
          arb_sub(term,term,tmp1,prec);
          neg_me=false;
        }
        else
        {
          arb_add(term,term,tmp1,prec);
          neg_me=true;
        }
      }
      arb_mul(term,term,a_pi_d[d],prec); // * (A pi)^d
      //printf("term(%" PRIu64 ",%" PRIu64 ") = ",n,d);arb_printd(term,10);printf("\n");
      exp_term(tmp,n,L,prec); // n/A = z, t0=0
      arb_mul(tmp1,term,tmp,prec); // (A Pi)^d * exp(Pi*r*n/(4*A)-Pi*n^2/A^2H^2)
      arb_mul(tmp,L->u_values_off[0][n*L->u_stride],tmp1,prec); // * Lambda(n/A)
      arb_add(res,res,tmp,prec);
      arb_mul(tmp,L->u_values_off[0][-n*L->u_stride],tmp1,prec); // * Lambda(-n/A)
      if(d&1)
        arb_neg(tmp,tmp);
      arb_add(res,res,tmp,prec);
    }
    arb_add_error(res,L->upsampling_error);
    return true;
  }

  Lerror_t do_rank(Lfunc *L)
  {
    if(!arb_contains_zero(L->u_values_off[0][0]))
    {
      arb_set(L->Lam_d,L->u_values_off[0][0]);
      if((L->rank==0)||(L->rank==DK))
      {
        L->rank=0;
        return ERR_SUCCESS;
      }
      else
        return ERR_CONFLICT_RANK;
    }
    int rank;
    for(rank=1;rank<=MAX_L;rank++)
    {
      df_zero(L->Lam_d,rank,L,L->wprec);
      if(!arb_contains_zero(L->Lam_d))
      {
        if((L->rank==rank)||(L->rank==DK))
        {
          L->rank=rank;
          return ERR_SUCCESS;
        }
        return ERR_CONFLICT_RANK;
      }
    }
    if(L->rank==DK) // can't leave it negative
      L->rank=0; // so guess at 0
    return ERR_NO_RANK;
  }


#ifdef __cplusplus
}
#endif
