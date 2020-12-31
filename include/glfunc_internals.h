#ifndef _GLFUNC_INTERNAL
#define _GLFUNC_INTERNAL
#define OUTPUT_RATIO (8) // we will analyse this portion of B
#define TURING_RATIO (16)
#define EXTRA_BITS (35) // extra bits of precision for convolves etc.
#define verbose (false)
#define BAD_64 (1LL<<62)

// how many integrals have I precomputed for Buthe's method?
#define MAX_MUI (100)
#define MAX_MUI_2 (200)
#define MAX_MU ((double) MAX_MUI)
#define MAX_MU_2 ((double) MAX_MUI_2)

#define MAX_L (10) // maximum differential allowed in upsampling

#include "inttypes.h"
#include "acb.h"
#include <stdbool.h>
#include "glfunc.h"

#ifdef __cplusplus
extern "C"{
#endif

  typedef struct{
    uint64_t degree;
    uint64_t conductor;
    double normalisation;
    double *mus;
    int64_t target_prec;
    arb_t zero_prec;
    arb_t zero_error;
    int64_t wprec; // working precision
    int64_t gprec; // precison used by g
    char *cache_dir;
    int self_dual;
    int rank;
    arb_t mu;
    arb_t nu;
    arb_t *nus;
    arb_t C;
    arb_t alpha;
    //arb_t k0; // no longer used
    //arb_t k2;
    double one_over_B;
    arb_t B;
    arb_t two_pi_by_B;
    arb_t pi;
    int64_t low_i;
    int64_t hi_i;
    uint64_t max_K;
    arb_t eq59;

    arb_t **Gs;

    // computation related
    uint64_t fft_N;
    uint64_t fft_NN;
    double A;
    arb_t arb_A;
    arf_t arf_A;
    arb_t one_over_A;
    arf_t arf_one_over_A;
    acb_t *G;
    acb_t *w; // twiddle factors for length fft_N
    acb_t *ww; // ditto for fft_NN
    arb_t *zeros[2];
    double eta;
    arb_t delta;
    arb_t exp_delta;
    acb_t *kres,*res,**skm;
    arb_t pre_ftwiddle_error;
    arb_t ftwiddle_error;

    arb_t buthe_Wf;
    arb_t buthe_Winf;
    arb_t buthe_Ws;
    arb_t buthe_b;
    arb_t buthe_sig1;
    arb_t buthe_C;
    arb_t buthe_h;
    arb_t buthe_ints[(MAX_R-1)*(2*MAX_MUI_2+1)];
    uint64_t buthe_M;


    arb_t one_over_root_N;
    arb_t sum_ans;
    acb_t epsilon;
    acb_t epsilon_sqr;
    acb_t *ans;
    uint64_t M;
    uint64_t M0;
    uint64_t allocated_M;
    double dc;

    bool nmax_called; // true if user/system has called Lfunc_nmax

    int64_t offset;

    uint64_t u_N;
    //arb_t *u_coshs;
    //arb_t *u_exps;
    arb_t u_H;
    arb_t u_pi_by_H2; //-pi/H^2
    arb_t u_A;
    arb_t u_one_over_A;
    arb_t *u_values[2];
    arb_t *u_values_off[2];
    uint64_t u_no_values;
    uint64_t u_no_values_off;
    uint64_t u_stride;
    arb_t u_pi_A;
    arb_t upsampling_error;

    arb_t Lam_d; // Lambda^(rank)(1/2)
    arb_t L_d; // L^(rank)(1/2)/rank!
  } Lfunc;

  // from glfunc_g.c
  Lerror_t compute_g(Lfunc *);

  // from acb_fft.c
  void acb_initfft(acb_t *w, uint64_t n, uint64_t prec);
  void acb_fft(acb_t *x, uint64_t n, acb_t *w, uint64_t prec);
  void acb_ifft(acb_t *x, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve1(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve2(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);

  // from error.c
  void abs_gamma(arb_t res, acb_t s, Lfunc *L, int64_t prec);
  void init_ftwiddle_error(Lfunc *L, int64_t prec);
  void complete_ftwiddle_error(Lfunc *L, int64_t prec);
  Lerror_t do_pre_iFFT_errors(Lfunc *L);
  bool M_error(arb_t res, arb_t x, Lfunc *L, int64_t prec);

  // from buthe.c
  void init_buthe(Lfunc *L, int64_t prec);
  void wf(Lfunc *L, uint64_t p, acb_poly_t fp1, acb_poly_t fp, int64_t prec);
  void buthe_Wf_error(Lfunc *L);
  Lerror_t buthe_check_RH(Lfunc *L);

  // from compute.c
  void lfunc_compute(Lfunc *L);

  //from upsample.c
  double upsample_error(long double M, long double H, long double h, long double A, double *mus, uint64_t r, uint64_t N, long double T, long double imz, uint64_t l);
  Lerror_t init_upsampling(Lfunc *L);
  bool newton(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec);
  bool upsample_stride(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec);
  Lerror_t arb_upsampling_error(arb_t res, double M,double H,double h,double A,double *mus,uint64_t r,uint64_t N,double T,arb_t imz, uint64_t l, arb_t pi, int64_t prec);

  //from zeros.c
  Lerror_t find_zeros(Lfunc *L, uint64_t side);

  // from rank.c
  //uint64_t guess_rank(Lfunc *L, uint64_t side, uint64_t prec);
  Lerror_t do_rank(Lfunc *L);

#ifdef __cplusplus
}
#endif
#endif
