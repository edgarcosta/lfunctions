#include "glfunc.h"
#include "glfunc_internals.h"

#undef verbose
#define verbose false
#ifdef __cplusplus
extern "C"{
#endif


#define POS (1)
#define NEG (2)
#define UNK (3)  
#define UP (1)
#define DOWN (2)
#define BOTH (3)
#define sign_t uint8_t
#define direction_t uint8_t

#define MAX_ZEROS (256)

  sign_t sign(arb_t x)
  {
    if(arb_contains_zero(x))
      return UNK;
    if(arb_is_positive(x))
      return POS;
    return NEG;
  }

  // which way is the curve going?
  direction_t direction(arb_t a, arb_t b, uint64_t prec)
  {
    static bool init=false;
    static arb_t tmp;
    if(!init)
    {
      init=true;
      arb_init(tmp);
    }
    arb_sub(tmp,a,b,prec);
    if(arb_contains_zero(tmp))
      return UNK;
    if(arb_is_negative(tmp)) // b>a
      return UP;
    return DOWN; // b<a
  }

  // binary chop
  bool isolate_zero(arb_t res, arb_t tt0, arb_t tt1, arb_t ff0, arb_t ff1, int8_t s0, int8_t s1 __attribute__((unused)), Lfunc *L, uint64_t side, uint64_t prec, uint64_t op_acc)
  {
    static bool init=false;
    static arb_t tmp1, tmp2, tmp3, t0, t1, f0, f1;
    if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(t0);
      arb_init(t1);
      arb_init(f0);
      arb_init(f1);

    }
    arb_set(t0,tt0);
    arb_set(t1,tt1);
    arb_set(f0,ff0);
    arb_set(f1,ff1);
    uint64_t iter=0;
    while(true)
    {
      arb_add(tmp1,t0,t1,prec);
      arb_mul_2exp_si(tmp1,tmp1,-1);

      // see if we are close enough yet
      arb_sub(tmp2,t1,t0,prec);
      arb_mul_2exp_si(tmp2,tmp2,op_acc+2); // +1 wasn't tight enough
      arb_sub_ui(tmp3,tmp2,1,prec);
      if(arb_is_negative(tmp3))
      {
        if(verbose) printf("Zero resolved in binary chop.\n");
        arb_set(res,tmp1);
        return true;
      }

      iter++;
      upsample_stride(tmp2,tmp1,L,side,prec);
      sign_t new_sign=sign(tmp2);
      if(new_sign==UNK)
      {
        if(verbose) printf("Indeterminate sign in upsampling on iteration %" PRIu64 ".\n",iter);
        arb_union(res,t0,t1,prec);
        if(verbose){printf("Best we could do was ");arb_printd(res,30);printf("\n");}
        return false; // flag to say not to full accuracy
      }
      if(new_sign!=s0) // change is between t0 and tmp1
      {
        arb_set(t1,tmp1);
        arb_set(f1,tmp2);
      }
      else // change is between tmp1 and t1
      {
        arb_set(t0,tmp1);
        arb_set(f0,tmp2);
      }
    }
  }

  // pass it an exact x=N*2^{-OP_ACC}, N an integer
  // check rigorously that there is a sign change between x\pm2^{-OP_ACC-1}
  // if so, then zero is within 2^{-OP_ACC-1} of x
  // this does Table Maker's dilemma
  bool confirm_zero(arb_t x, Lfunc *L, int64_t prec, arb_t zero_prec, uint64_t side)
  {
    static bool init=false;
    static arb_t fleft,fright,tmp;
    if(!init)
    {
      init=true;
      arb_init(fleft);
      arb_init(fright);
      arb_init(tmp);
    }
    if(verbose){printf("Confirming zero at ");arb_printd(x,100);printf("\n");}
    arb_sub(tmp,x,zero_prec,prec);
    if(arb_is_negative(tmp))
    {
      side^=1;
      arb_neg(tmp,tmp);
      //printf("computing F(");arb_printd(tmp,10);printf(") on side %" PRIu64 "\n",side);
      if(!upsample_stride(fleft, tmp, L, side, prec))
        return(false);
      //printf("Left = ");arb_printd(fleft,20);printf("\n");
      side^=1;
    }
    else
    {
      //printf("computing F(");arb_printd(tmp,10);printf(") on side %" PRIu64 "\n",side);
      if(!upsample_stride(fleft, tmp, L, side, prec))
        return(false);
      //printf("Left = ");arb_printd(fleft,20);printf("\n");
    }
    arb_add(tmp,x,zero_prec,prec);
    //printf("computing F(");arb_printd(tmp,10);printf(") on side %" PRIu64 "\n",side);
    if(!upsample_stride(fright, tmp, L, side, prec))
      return(false);
    //printf("Right = ");arb_printd(fright,20);printf("\n");
    if((sign(fleft)&sign(fright))==0)
      return true;
    if(verbose)
    {
      printf("Problem in confirm zero.\nf(t-delta)=");
      arb_printd(fleft,20);
      printf("\nf(t+delta)=");
      arb_printd(fright,20);
      printf("\n");
    }
    return false;
  }

  // snap rho to N*2^{-OP_ACC} where N is an integer
  bool snap_point(arb_t res, arb_t rho, uint64_t prec, uint64_t op_acc)
  {
    static bool init=false;
    static arb_t tmp;
    static mpfr_t a,b;
    static mpz_t z;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      mpfr_init2(a,prec);
      mpfr_init2(b,prec);
      mpz_init(z);
    }
    arb_mul_2exp_si(tmp,rho,op_acc);
    arb_get_interval_mpfr(a,b,tmp);
    mpfr_add(a,a,b,MPFR_RNDN);
    mpfr_div_ui(b,a,2,MPFR_RNDN);
    mpfr_get_z(z,b,MPFR_RNDN);
    mpfr_set_z(a,z,MPFR_RNDN);
    arb_set_interval_mpfr(tmp,a,a,prec);
    arb_mul_2exp_si(res,tmp,-op_acc);
    return true;
  }

  // called with suspected stationary point between m-1, m and m+1
  // if !isolate_p, just confirm the two zeros to minimal precision
  error_t stat_point(arb_t z1, arb_t z2, uint64_t m, Lfunc *L, uint64_t side, uint64_t prec, bool isolate_p)
  {
    error_t ecode=ERR_SUCCESS;
    static bool init=false;
    static arb_t t0,t1,f0,f1,t2,f2,t01,f01,t12,f12,tmp;
    if(!init)
    {
      init=true;
      arb_init(t0);
      arb_init(t1);
      arb_init(t2);
      arb_init(f0);
      arb_init(f1);
      arb_init(f2);
      arb_init(t01);
      arb_init(t12);
      arb_init(f01);
      arb_init(f12);
      arb_init(tmp);
    }

    arb_mul_ui(t0,L->one_over_A,m-1,prec);
    arb_mul_ui(t1,L->one_over_A,m,prec);
    arb_mul_ui(t2,L->one_over_A,m+1,prec);
    arb_mul_ui(tmp,t0,L->degree,prec);
    arb_mul(tmp,tmp,L->pi,prec);
    arb_mul_2exp_si(tmp,tmp,-2);
    arb_exp(tmp,tmp,prec);
    arb_mul(f0,L->u_values_off[side][m-1],tmp,prec);
    arb_mul_ui(tmp,t1,L->degree,prec);
    arb_mul(tmp,tmp,L->pi,prec);
    arb_mul_2exp_si(tmp,tmp,-2);
    arb_exp(tmp,tmp,prec);
    arb_mul(f1,L->u_values_off[side][m],tmp,prec);
    arb_mul_ui(tmp,t2,L->degree,prec);
    arb_mul(tmp,tmp,L->pi,prec);
    arb_mul_2exp_si(tmp,tmp,-2);
    arb_exp(tmp,tmp,prec);
    arb_mul(f2,L->u_values_off[side][m+1],tmp,prec);
    sign_t s=sign(f0);

    while(true)
    {
      //if(verbose){arb_printd(t0,20);printf(" ");arb_printd(t1,20);printf(" ");arb_printd(t2,20);printf("\n");}
      //if(verbose){arb_printd(f0,20);printf(" ");arb_printd(f1,20);printf(" ");arb_printd(f2,20);printf("\n");}
      arb_add(t01,t0,t1,prec);
      arb_mul_2exp_si(t01,t01,-1);
      if(!upsample_stride(f01,t01,L,side,prec))
        return ecode|ERR_STAT_POINT;
      //if(verbose){printf("t01 = ");arb_printd(t01,20);printf("\n");}
      //if(verbose){printf("f01 = ");arb_printd(f01,20);printf("\n");}
      sign_t s01=sign(f01);
      if(s01==UNK)
        return ecode|ERR_DBL_ZERO;
      if(s01!=s)
      {
        if(!isolate_p)
        {
          arb_union(z1,t0,t01,prec);
          arb_union(z2,t01,t1,prec);
          return ecode;
        }
        if(!isolate_zero(tmp,t0,t01,f0,f01,s,s01,L,side,prec,L->target_prec))
        {
          ecode|=ERR_ZERO_PREC;
          arb_set(z1,tmp);
        }
        else
        {
          if(!snap_point(z1,tmp,prec,L->target_prec))
            return ecode|ERR_STAT_POINT;
          if(!confirm_zero(z1,L,prec,L->zero_prec,side))
            return ecode|ERR_STAT_POINT;
        }
        if(!isolate_zero(tmp,t01,t1,f01,f1,s01,s,L,side,prec,L->target_prec))
        {
          ecode|=ERR_ZERO_PREC;
          arb_set(z2,tmp);
        }
        else
        {
          if(!snap_point(z2,tmp,prec,L->target_prec))
            return ecode|ERR_STAT_POINT;
          if(!confirm_zero(z2,L,prec,L->zero_prec,side))
            return ecode|ERR_STAT_POINT;
        }
        return ecode;
      }
      direction_t left=direction(f0,f01,prec);
      direction_t right=direction(f01,f1,prec);
      if((left==UNK)||(right==UNK))
        return ecode|ERR_DBL_ZERO;
      if(left!=right)
      {
        arb_set(t2,t1);
        arb_set(f2,f1);
        arb_set(t1,t01);
        arb_set(f1,f01);
        continue;
      }

      arb_add(t12,t1,t2,prec);
      arb_mul_2exp_si(t12,t12,-1);
      upsample_stride(f12,t12,L,side,prec);
      //printf("right middle = ");arb_printd(f12,20);printf("\n");
      sign_t s12=sign(f12);
      if(s12==UNK)
        return ecode|ERR_DBL_ZERO;
      if(s12!=s)
      {
        if(!isolate_p)
        {
          arb_union(z1,t1,t12,prec);
          arb_union(z2,t12,t2,prec);
          return ecode;
        }
        if(!isolate_zero(tmp,t1,t12,f1,f12,s,s12,L,side,prec,L->target_prec))
        {
          ecode|=ERR_ZERO_PREC;
          arb_set(z1,tmp);
        }
        else
        {
          if(!snap_point(z1,tmp,prec,L->target_prec))
            return ecode|ERR_STAT_POINT;
          if(!confirm_zero(z1,L,prec,L->zero_prec,side))
            return ecode|ERR_STAT_POINT;
        }
        if(!isolate_zero(tmp,t12,t2,f12,f2,s12,s,L,side,prec,L->target_prec))
        {
          ecode|=ERR_ZERO_PREC;
          arb_set(z2,tmp);
        }
        else
        {
          if(!snap_point(z2,tmp,prec,L->target_prec))
            return ecode|ERR_STAT_POINT;
          if(!confirm_zero(z2,L,prec,L->zero_prec,side))
            return ecode|ERR_STAT_POINT;
        }
        return ecode;
      }
      left=direction(f1,f12,prec);
      right=direction(f12,f2,prec);
      if((left==UNK)||(right==UNK))
        return ecode|ERR_DBL_ZERO;
      if(left!=right)
      {
        arb_set(t0,t1);
        arb_set(f0,f1);
        arb_set(t1,t12);
        arb_set(f1,f12);
        continue;
      }
      else
      {
        arb_set(t0,t01);
        arb_set(f0,f01);
        arb_set(t2,t12);
        arb_set(f2,f12);
        continue;
      }
    }
    return ecode;
  }

  error_t find_zeros(Lfunc *L, uint64_t side)
  {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,tmp3,t0,z1,z2;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(t0);
      arb_init(z1);
      arb_init(z2);
    }
    int64_t prec=L->wprec;
    bool stat_points=true;
    for(uint64_t z=0;z<MAX_ZEROS;z++)
      arb_zero(L->zeros[side][z]);
    uint64_t op_acc=L->target_prec;
    uint64_t count=0;

    uint64_t n=0;
    sign_t this_sign,last_sign=sign(L->u_values_off[side][0]);
    direction_t this_dir,last_dir;
    error_t ecode=ERR_SUCCESS;
    last_dir=direction(L->u_values_off[side][0],L->u_values_off[side][1],prec);
    if(last_dir==UNK)
    {
      printf("Unknown direction at start of data.\n");
      return ecode|ERR_NO_DATA;
    }

    if(last_sign==UNK) // central zero(s)
    {
      n++; // skip the central pt.
      last_sign=sign(L->u_values_off[side][n]);
      if(last_sign==UNK) // that was your last chance
        return ecode|ERR_NO_DATA;
    }

    // now start searching for zeros and stat pts.
    while(true)
    {
      n++;
      if(n>L->fft_NN/OUTPUT_RATIO+L->fft_NN/TURING_RATIO)
        return ecode;
      this_sign=sign(L->u_values_off[side][n]);
      if(this_sign==UNK)
      {
        //printf("Indeterminate sign.\n"); // time to give up
        return ecode|ERR_SOME_DATA;
      }
      if(this_sign!=last_sign) // found a zero between n and n-1
      {
        arb_mul_ui(tmp1,L->one_over_A,n-1,prec);
        arb_mul_ui(tmp2,L->one_over_A,n,prec);
        if(verbose)
        {
          printf("zero found between ");arb_printd(tmp1,20);
          printf(" and ");arb_printd(tmp2,20);printf("\n");
        }
        if(n<L->fft_NN/OUTPUT_RATIO)// need to isolate it properly
          // always isolate to reduce chances of straddling start of turing zone
        {
          arb_add(tmp3,tmp1,tmp2,prec);
          arb_mul_2exp_si(tmp3,tmp3,-1);
          if(!newton(tmp,tmp3,L,side,prec)||!snap_point(tmp3,tmp,prec,op_acc)||!confirm_zero(tmp3,L,prec,L->zero_prec,side))
          {
            if(verbose) printf("Resorting to binary chop.\n");
            if(!isolate_zero(tmp,tmp1,tmp2,L->u_values_off[side][n-1],L->u_values_off[side][n],last_sign,this_sign,L,side,prec,op_acc))
            {
              if(verbose){printf("Setting zeros number %" PRIu64 " to ",count);arb_printd(tmp,20);printf("\n");}
              arb_set(L->zeros[side][count++],tmp);
              ecode|=ERR_ZERO_PREC;
            }
            else  // isolate zero got to +/- target_prec
            {
              if(!snap_point(tmp3,tmp,prec,op_acc)) return ecode|=ERR_ZERO_ERROR;
              if(!confirm_zero(tmp3,L,prec,L->zero_prec,side))
              {
                printf("Failed to confirm zero near ");
                arb_printd(tmp3,100);
                printf("\n");
                return ecode|=ERR_ZERO_ERROR;
              }
              if(verbose){printf("Setting zeros number %" PRIu64 " to ",count);arb_printd(tmp3,20);printf("\n");}
              arb_set(L->zeros[side][count],tmp3);
              //arb_add_error(L->zeros[side][count++],L->zero_prec);
              arb_add_error_2exp_si(L->zeros[side][count++],-L->target_prec-1);
            }
          }
          else // newton iteration got us to +/- target_prec
          {
            if(verbose){printf("setting zero number %" PRIu64 " to ",count);arb_printd(tmp3,200);printf("\n");}
            arb_set(L->zeros[side][count],tmp3);
            //arb_add_error(L->zeros[side][count++],L->zero_prec);
            arb_add_error_2exp_si(L->zeros[side][count++],-L->target_prec-1);
          }
        }
        else // don't isolate it any more, just record it
        {
          arb_mul_ui(tmp1,L->one_over_A,n-1,prec);
          arb_mul_ui(tmp2,L->one_over_A,n,prec);
          arb_union(L->zeros[side][count],tmp1,tmp2,prec);
          if(verbose){printf("setting zero number %" PRIu64 " to ",count);arb_printd(L->zeros[side][count],20);printf("\n");}
          count++;
        }
        last_sign=this_sign;
        last_dir=direction(L->u_values_off[side][n-1],L->u_values_off[side][n],prec);
        if(last_dir==UNK) // unkown dirs in region. give up on stat points
          stat_points=false;
        continue;
      }

      if(!stat_points) // given up on looking for stationary points
        continue;
      this_dir=direction(L->u_values_off[side][n-1],L->u_values_off[side][n],prec);
      if(this_dir==UNK)
      {
        stat_points=false; // time to give up on stat points
        continue;
        //printf("Unknown direction in find_zeros.\n");
        //return 0;
      }

      if((this_dir&last_dir)==0) // change in direction
      {
        if(((last_dir==UP)&&(this_sign==NEG))||((last_dir==DOWN)&&(this_sign==POS)))
        {
          if(verbose) printf("Stationary point detected.\n");
          ecode|=stat_point(z1,z2,n-1,L,side,prec,n<=L->fft_NN/OUTPUT_RATIO);
          if(fatal_error(ecode))
            return ecode;
          if(verbose)
          {
            printf("Stat zero found at ");arb_printd(z1,20);
            printf("\nand                ");arb_printd(z2,20);
            printf("\n");fflush(stdout);
          }
          if(verbose){printf("setting zeros %" PRIu64 " and %" PRIu64 " to ",count,count+1);
            arb_printd(z1,20);printf(" ");arb_printd(z2,20);printf("\n");}
          arb_set(L->zeros[side][count],z1);
          arb_set(L->zeros[side][count+1],z2);
          count+=2;
        }

        last_dir=this_dir;
      }
      last_dir=this_dir;
    }
  }

#ifdef __cplusplus
}
#endif
