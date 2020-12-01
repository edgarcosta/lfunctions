// Copyright Edgar Costa 2019
// See LICENSE file for license details.
//
/*
 * Computes the L-function data for a smalljac curve
 *
 *
 * Usage:
 *
 * Input file:
 * label:cond:curve_str:bad_euler_factors
 * bad_euler_factors = [[Norm(p1), Norm(p2), ...], [Lp1, Lp2, ..]]
 * where pi is a prime that divides the discriminant of the curve
 * For example:
11.a:11:[0, -1, 1, 0, 0]:[[11],[[1,-1]]]
109.a:109:[1, -1, 0, -8, -7]:[[109],[[1,-1]]]
5077.a:5077:[0, 0, 1, -7, 6]:[[5077],[[1,1]]]
112.c:112:[0,-1,0,-43688,3529328]:[[2,7],[[1],[1,1]]]
7406.a:7406:[1,0,1,-276,-3586]:[[2,7,23],[[1,1],[1,1],[1]]]
277.a:277:[-x^2-x, x^3+x^2+x+1]:[[277], [[1, -7, 269, 277]]]
25913.a.25913.1:25913:[x^3 - x^2 - 2*x, x^3 + x + 1]:[[25913], [[1, -62, 25974, -25913]]]
and elliptic curves over quadratic fields, note that one now passes thed bad factors via the prime norm:
76.1-a:1900:[a,0,a,a,0]/(a^2-a-1):[[4,19],[[1,0,1],[1,1]]]
2.2.92.1-98.1-f:7406:[a,0,a,-1,0]/(a^2-23):[[2,7,7,23],[[1,1],[1,1],[1,1],[1,0,23]]]
2.0.3.1-120000.1.b:1080000:[0,-1,0,7,-3]/(a^2-a+1):[[4,3,25],[[1],[1,1],[1]]]

 *
 * Output file:
 * label:root number:rank:leading term taylor:10 zeros:plot delta:plot yvalues
 */

#ifdef NOSMALLJAC

#include <stdio.h>
int main() {
  printf("Binary not compiled, no smalljac or ff_poly detected at compilation time");
  return 0;
}

#else

#define __STDC_FORMAT_MACROS
#include <chrono>
#include <cstdint>
#include <cstring>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpzxx.h>
#include <flint/fmpz_polyxx.h>
#include <flint/nmod_poly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/ulong_extras.h>
#include <flint/perm.h>
#include <primesieve.h>
#include <smalljac.h>

#include "glfunc.h"
#include "examples_tools.h"

using std::strcpy;
using flint::fmpz_polyxx;

typedef std::chrono::time_point<std::chrono::system_clock> SystemTime ;


typedef struct {
  // some string
  string label;

  // the conductor of the jacobian
  int64_t conductor;


  // the string input for smalljac
  string sj_str;
  smalljac_curve_t sj_curve;
  int sj_err;
  int genus;
  int nf_degree;

  // the bad local factors
  multimap<int64_t, vector<fmpzxx>> bad_factors;


  // nf field stuff
  fmpz_polyxx nf_poly;
  fmpzxx nf_disc;
  unsigned long lastq;
  fmpz_polyxx last_local_factor;
  vector<bool> touched; // the euler factors already used



  // L-function
  Lerror_t ecode;
  Lfunc_t L;
  double* mus;
  int degree;
} curve;



istream &operator>>(istream &is, curve &o)
{
  for(size_t i = 0; i < 4; ++i){
    string s;
    if(not getline(is, s, ':'))
      throw_line("missing field"s);
    stringstream ss(s);
    char* cstr;
    std::string nf_poly_string;
    std::stringstream nf_poly_string_ss;
    switch(i) {
      case 0:
        if(!(ss >> o.label))
          throw_line("bad label"s);
        break;
      case 1:
        if(!(ss >> o.conductor))
          throw_line("bad conductor"s);
        break;
      case 2:
        o.sj_str = ss.str();
        cstr = new char[o.sj_str.length()+1];
        strcpy(cstr, o.sj_str.c_str());
        o.sj_err = 0;
        o.sj_curve = smalljac_curve_init(cstr, &o.sj_err);
        delete[] cstr;
        if(o.sj_err != 0) {
          stringstream message;
          message << "smalljac error = "<< o.sj_err << " while processing " << o.sj_str;
          throw_line(message.str());
        }
        o.genus = smalljac_curve_genus(o.sj_curve);
        o.nf_degree = smalljac_curve_nf_degree(o.sj_curve);
        // smalljac only handles quadratic number fields for elliptic curves
        assert_print(o.nf_degree, <=, 2);
        if(o.nf_degree == 2) {
          nf_poly_string = std::string(smalljac_curve_nf(o.sj_curve));
          // remove ( )
          nf_poly_string.erase(0, 1);
          nf_poly_string.erase(nf_poly_string.size() - 1, 1);
          nf_poly_string_ss = std::stringstream(nf_poly_string);
          nf_poly_string_ss >> o.nf_poly;
          fmpz_poly_discriminant(o.nf_disc._fmpz(), o.nf_poly._poly());
          // make it maximal at 2
          while( o.nf_disc.divisible(4) ) { o.nf_disc = o.nf_disc.divexact(4); }
          if ( (o.nf_disc % fmpzxx(4)) != 1 ) o.nf_disc *= 4;
        }
        o.degree = 2*o.genus*o.nf_degree;
        o.mus = new double[o.degree];
        o.lastq = 1;
        for(int i = 0; i < o.genus*o.nf_degree; ++i) {
          o.mus[i] = 0;
          o.mus[i + o.genus*o.nf_degree] = 1;
        }

        o.L = Lfunc_init(o.degree, uint64_t(o.conductor), 0.5, o.mus, &o.ecode);
        break;
      case 3:
        if(!(ss >> o.bad_factors))
          throw_line("bad input for bad local factors"s);
        break;
      default:
        throw_line("too many fields in the line!"s);
    }
  }
  return is;
}

void curve_clear(curve &C) {
  smalljac_curve_clear(C.sj_curve);
  delete[] C.mus;
  Lfunc_clear(C.L);
}

int smalljac_callback(
     __attribute__((unused)) smalljac_curve_t c,
    unsigned long q,		// prime (or prime power) in [start,end]
    int good,						// 1 if good reduction, 0 if bad
    long a[],						// n coefficients of L_q(T) with a[0] = a_1
    int n,							// either 1 or g for good p, 0 for bad p
    void *arg){  				// forwarded arg from caller

	curve *C = (curve *)arg;
  acb_poly_t local_factor;
  acb_poly_init(local_factor);
  fmpz_polyxx local_factor_zz;

  if( good ) {
    local_factor_zz.fit_length(n == C->genus ? 2*C->genus + 1 : 1 + n);
    local_factor_zz.set_coeff(0, 1);
    for(int i = 0; i < n; ++i)
      local_factor_zz.set_coeff(i + 1, a[i]);
    // complete with the functional equation
    if( n == C->genus) {
      fmpzxx qn(q);
      for(int i = C->genus + 1; i <= 2*C->genus; ++i) {
        local_factor_zz.set_coeff(i, local_factor_zz.get_coeff(2*C->genus - i)*qn);
        qn *= q;
      }
    }
  } else {
    auto it = C->bad_factors.find(int64_t(q));
    if(it == C->bad_factors.end()) {
      stringstream message;
      message << "local factor for q = "<< q << " not found!";
      throw_line(message.str());
    }
    local_factor_zz.fit_length(it->second.size());
    for(int i = 0; i < it->second.size(); ++i)
      local_factor_zz.set_coeff(i, it->second[i]);
    C->bad_factors.erase(it); // we will no longer use it
  }
  bool use_lpoly = true;
  if(C->nf_degree == 2){
    if(n_is_square(q)) {
        q = n_sqrt(q); // so we call Lfunc_use_lpoly accordingly
        fmpz_polyxx local_factor_zz2(2*local_factor_zz.degree() + 1);
        for(size_t i=0; i <= 2*local_factor_zz.degree(); ++i) {
          if(i%2 == 0)
            local_factor_zz2.set_coeff(i, local_factor_zz.get_coeff(i/2));
        }
        local_factor_zz = local_factor_zz2;
    } else { // p is split, as smalljac doesn't handle ramified primes in the monic order
      if(C->lastq == q) {
        local_factor_zz *= C->last_local_factor;
      } else {
        C->last_local_factor = local_factor_zz;
        use_lpoly = false;
      }
    }
    C->lastq = q;
  }
  if(use_lpoly) {
    acb_poly_fit_length(local_factor, local_factor_zz.degree() + 1);
    _acb_poly_set_length(local_factor, local_factor_zz.degree() + 1);
    for(long i = 0; i <= local_factor_zz.degree(); ++i) {
      // acb_poly_get_coeff_ptr(local_factor, i) = local_factor->coeffs + i
      acb_set_fmpz(local_factor->coeffs + i, local_factor_zz.coeff(i)._fmpz());
    }
    Lfunc_use_lpoly(C->L, q, local_factor);
    C->touched[q] = true;
  }
  acb_poly_clear(local_factor);
  return true;
}

// what does this return?
long populate_local_factors(curve &C) {
  C.touched = vector<bool>(Lfunc_nmax(C.L) + 1, false);
  long res = smalljac_Lpolys(C.sj_curve, 1, Lfunc_nmax(C.L), 0, smalljac_callback, &C);
  if(C.nf_degree == 2) { // deal with ramified primes
    acb_poly_t local_factor;
    acb_poly_init(local_factor);
    for(const auto &it : C.bad_factors){
      assert_print(C.nf_disc % fmpzxx(it.first), ==, 0);
      if( it.first <= (long) Lfunc_nmax(C.L) ) {
        acb_poly_fit_length(local_factor, it.second.size());
        _acb_poly_set_length(local_factor, it.second.size());
        for(long i = 0; i < it.second.size(); ++i) {
          // acb_poly_get_coeff_ptr(local_factor, i) = local_factor->coeffs + i
          acb_set_fmpz(local_factor->coeffs + i, it.second[i]._fmpz());
        }
        Lfunc_use_lpoly(C.L, it.first, local_factor);
        C.touched[it.first] = true;
      }
    }
    acb_poly_one(local_factor);
    primesieve_iterator it;
    primesieve_init(&it);
    uint64_t p=0;
    while((p=primesieve_next_prime(&it)) <= Lfunc_nmax(C.L)){
      if( not C.touched[p] ) Lfunc_use_lpoly(C.L, p, local_factor);
    }
    primesieve_free_iterator(&it);
    acb_poly_clear(local_factor);
  }
  return res;
}



ostream& operator<<(ostream &s, curve &C) {
  Lfunc_t &L = C.L;

  s << C.label <<":";
  // root number
  s << Lfunc_epsilon(L) <<":";
  // r = rank
  s << Lfunc_rank(L) << ":";
  // L(1/2)^r / r! as arb
  s << Lfunc_Taylor(L) << ":";
  // first zeros
  arb_srcptr zeros = Lfunc_zeros(L, 0);
  s << "[";
  for(size_t i = 0; i < 10; ++i) {
    s << zeros + i;
    if( i < 9 )
      s << ", ";
    else
      s << "]";
  }
  s << ":";
  Lplot_t *Lpp = Lfunc_plot_data(L, 0, 10.0, 500);
  s << Lpp;
  Lfunc_clear_plot(Lpp);
  return s;
}




int main (int argc, char**argv)
{
  try {
    assert_print(argc, ==, 3);
    printf("Input: %s\n", argv[1]);
    printf("Output: %s\n", argv[2]);

    ifstream input(argv[1]);
    ofstream output(argv[2]);
    string   line;

    int r = 0;

    while(std::getline(input, line)) {
      SystemTime start(std::chrono::system_clock::now());
      std::time_t startt = std::chrono::system_clock::to_time_t(start);
      cout << "Date:   \t" <<  std::put_time(std::localtime(&startt), "%F %T") << endl;

      curve C;
      Lerror_t &ecode = C.ecode;
      Lfunc_t &L = C.L;


      // read a line
      stringstream linestream(line);
      linestream >> C;
      cout << "Starting:\t"<<C.label<<endl;

      // we need all the local factors p <= target_M
      uint64_t target_M = Lfunc_nmax(L);

      cout <<"\tusing p <= " << target_M << endl;

      // populate local factors
      populate_local_factors(C);

      // do the computation
      ecode |= Lfunc_compute(L);
      if(fatal_error(ecode)) {
        fprint_errors(stderr, ecode);
        std::abort();
      }

      printf("\tRank = %" PRIu64 "\n",Lfunc_rank(L));
      printf("\tEpsilon = ");acb_printd(Lfunc_epsilon(L),20);printf("\n");
      printf("\tFirst non-zero Taylor coeff = ");arb_printd(Lfunc_Taylor(L), 20);printf("\n");
      printf("\tFirst zero = ");arb_printd(Lfunc_zeros(L, 0), 20);printf("\n");


      output << C << endl;
      // print any warnings collected along the way
      if( ecode != ERR_SUCCESS )
      {
        cerr << "Begin warnings for " << C.label << endl;
        fprint_errors(stderr, ecode);
        cerr << "End warnings for " << C.label << endl;
        r++;
      }

      SystemTime end(std::chrono::system_clock::now());
      std::time_t endt = std::chrono::system_clock::to_time_t(end);
      double walltime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      cout << "Done:   \t"<< C.label << "\ttook ";
      cout << std::setw(6) << std::setfill(' ')  << std::fixed << std::setprecision(2) << walltime/1000 << "s"<< endl;
      cout << "Date:   \t" <<  std::put_time(std::localtime(&endt), "%F %T") << endl << endl;

      //free memory
      curve_clear(C);
    }
    flint_cleanup();
    return r;
  } catch( const std::exception & ex ) {
     cerr << "Uncaught exception: " <<ex.what() << endl;
     std::abort();
  }
}

#endif

