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
277.a:277:[-x^2-x, x^3+x^2+x+1]:[[277], [[1, -7, 269, 277]]]
25913.a.25913.1:25913:[x^3 - x^2 - 2*x, x^3 + x + 1]:[[25913], [[1, -62, 25974, -25913]]]
and in the future:
76.1-a:1900:[a,0,a,a,0]/(a^2-a-1):[[4,19],[[1,0,1],[1,1]]]

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
#include <flint/fmpzxx.h>
#include <flint/nmod_poly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/ulong_extras.h>
#include <flint/perm.h>
#include <smalljac.h>

#include "glfunc.h"
#include "examples_tools.h"

using std::strcpy;

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


  // L-function
  Lerror_t ecode;
  Lfunc_t L;
  double* mus;
  int degree;
} curve;


istream & operator>>(istream &is, curve &o)
{
  for(size_t i = 0; i < 4; ++i){
    string s;
    if(not getline(is, s, ':'))
      throw_line("missing field"s);
    stringstream ss(s);
    char* cstr;

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
        // for now no number fields
        assert_print(o.nf_degree, ==, 1);
        o.degree = 2*o.genus*o.nf_degree;
        o.mus = new double[o.degree];
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
  vector<fmpzxx> local_factor_zz;
  if( good ) {
    local_factor_zz.resize(n == C->genus ? 2*C->genus + 1 : 1 + n);
    local_factor_zz[0] = 1;
    for(int i = 0; i < n; ++i)
      local_factor_zz[i + 1] = a[i];
    // complete with the functional equation
    if( n == C->genus) {
      fmpzxx qn(q);
      for(int i = C->genus + 1; i <= 2*C->genus; ++i) {
        local_factor_zz[i] = local_factor_zz[2*C->genus - i]*qn;
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
    local_factor_zz = it->second;
    C->bad_factors.erase(it); // we will no longer use it
  }
  //if( C->nf_degree > 1 and n_is_square(q)) {

  //  q = n_sqrt(q);
  //  vector<fmpzxx> local_factor_zz2(2*local_factor_zz.size() - 1);
  //  for(size_t i=0; i < local_factor_zz2.size(); ++i) {
  //    if(i%2) { // == 1
  //      local_factor_zz2[i] = 0;
  //    } else {
  //      local_factor_zz2[i] = local_factor_zz[i/2];
  //    }
  //  }
  //  local_factor_zz.swap(local_factor_zz2);
  //}
  acb_poly_fit_length(local_factor, local_factor_zz.size());
  _acb_poly_set_length(local_factor, local_factor_zz.size());
  for(size_t i = 0; i < local_factor_zz.size(); ++i) {
    // acb_poly_get_coeff_ptr(local_factor, i) = local_factor->coeffs + i
    acb_set_fmpz(local_factor->coeffs + i, local_factor_zz[i]._fmpz());
  }
  Lfunc_use_lpoly(C->L, q, local_factor);
  acb_poly_clear(local_factor);
  return true;
}

// what does this return?
long populate_local_factors(curve &C) {
  print(Lfunc_nmax(C.L));
  return smalljac_Lpolys(C.sj_curve, 1, Lfunc_nmax(C.L), 0, smalljac_callback, &C);
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

      cout <<"using p <= " << target_M << endl;

      // populate local factors
      populate_local_factors(C);

      // do the computation
      ecode |= Lfunc_compute(L);
      if(fatal_error(ecode)) {
        fprint_errors(stderr, ecode);
        std::abort();
      }

      printf("Rank = %" PRIu64 "\n",Lfunc_rank(L));
      printf("Epsilon = ");acb_printd(Lfunc_epsilon(L),20);printf("\n");
      printf("First non-zero Taylor coeff = ");arb_printd(Lfunc_Taylor(L), 20);printf("\n");


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
      cout << "Date:   \t" <<  std::put_time(std::localtime(&endt), "%F %T") << endl;
      cout << "Done:   \t"<< C.label << "\ttook ";
      cout << std::setw(6) << std::setfill(' ')  << std::fixed << std::setprecision(2) << walltime/1000 << "s"<< endl << endl;

      //free memory
      curve_clear(C);
    }
    return r;
  } catch( const std::exception & ex ) {
     cerr << "Uncaught exception: " <<ex.what() << endl;
     std::abort();
  }
}

#endif

