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
11.a:1:11:[0, -1, 1, 0, 0]:[[11],[[1,-1]]]
11.a:2:121:[0, -1, 1, 0, 0]:[[11],[[1,-1]]]
11.a:3:1331:[0, -1, 1, 0, 0]:[[11],[[1,-1]]]
109.a:1:109:[1, -1, 0, -8, -7]:[[109],[[1,-1]]]
5077.a:1:5077:[0, 0, 1, -7, 6]:[[5077],[[1,1]]]
112.c:1:112:[0,-1,0,-43688,3529328]:[[2,7],[[1],[1,1]]]
7406.a:1:7406:[1,0,1,-276,-3586]:[[2,7,23],[[1,1],[1,1],[1]]]
277.a:1:277:[-x^2-x, x^3+x^2+x+1]:[[277], [[1, -7, 269, 277]]]
25913.a.25913.1:1:25913:[x^3 - x^2 - 2*x, x^3 + x + 1]:[[25913], [[1, -62, 25974, -25913]]]
and elliptic curves over quadratic fields, note that one now passes the bad factors via the prime norm:
2.2.5.1-76.1-a:1:1900:[a,0,a,a,0]/(a^2-a-1):[[4,5,19],[[1,1],[1,-1,5],[1,1]]]
2.2.92.1-98.1-f:1:829472:[a,0,a,-1,0]/(a^2-23):[[2,7,7,23],[[1,1],[1,1],[1,1],[1,0,23]]]
2.0.3.1-120000.1.b:1:1080000:[0,-1,0,7,-3]/(a^2-a+1):[[4,3,25],[[1],[1,1],[1]]]

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
#define special_values_size 2 // implies computing L(1) ... L(special_values_size)
#include <array>
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
#include "glfunc_internals.h" // so we can set self-dual = yes
#include "tracehash.h"
#include "examples_tools.h"

using std::array;
using std::strcpy;
using flint::fmpzxx;
using flint::fmpqxx;
using flint::fmpz_polyxx;

#define M61UI		0x1FFFFFFFFFFFFFFFUL		// this should raise a compile-time error on 32-bit machines


// i -> p_i
int16_t* primes;
// p_i -> i
map<int16_t, int16_t> primes_dict;

typedef std::chrono::time_point<std::chrono::system_clock> SystemTime;

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

  // Sym^d(L(C))
  int symdegree;

  //  k -> L(k+1)
  array<acb_struct, special_values_size> special_values;


  // nf field stuff
  fmpz_polyxx nf_poly;
  fmpzxx nf_disc;
  unsigned long lastq;
  fmpz_polyxx last_local_factor;
  vector<bool> touched; // the euler factors already used

  //stores ap, p < 2^13
  array<fmpzxx, APHASH_MAXPI> ap;

  //2nd moment approximation
  double second_moment;

  //trace hash
  long tracehash;


  // L-function
  Lerror_t ecode;
  Lfunc_t L;
  double* mus;
  int degree;
} curve;

istream &operator>>(istream &is, curve &o)
{
  // initialize ap
  for(auto &elt: o.ap)
    elt = 0;
  // initialize sum for second moment
  o.second_moment = 0;
  for(size_t i = 0; i < 5; ++i){
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
        if(!(ss >> o.symdegree))
          throw_line("bad sym degree"s);
        break;
      case 2:
        if(!(ss >> o.conductor))
          throw_line("bad conductor"s);
        break;
      case 3:
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
        if( o.genus > 1 or o.nf_degree > 1) {
          // we have only hard coded the euler factor formulas for symmetric powers for elliptic curves over Q
          assert_print(o.symdegree, ==, 1);
        }
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
        if(o.symdegree > 1) {
          //slong newdeg = flint::bin(ulong(o.degree + o.symdegree - 1), ulong(o.symdegree)).to<slong>();
          //o.degree = newdeg;
          o.degree = o.symdegree + 1;
          assert_print(o.degree, <=, MAX_DEGREE);
        }
        o.mus = new double[o.degree];
        o.lastq = 1;
        if(o.symdegree == 1) {
          for(int i = 0; i < o.genus*o.nf_degree; ++i) {
            o.mus[i] = 0;
            o.mus[i + o.genus*o.nf_degree] = 1;
          }
        } else {
          assert_print(o.genus, ==, 1);
          int u = ceil(o.symdegree * 0.5);
          for(int i = 0; i < u; ++i) {
            o.mus[i] = -i;
            o.mus[i + u] = -i + 1;
          }
          if((o.symdegree % 2) == 0)
            o.mus[o.symdegree] = -2*floor(u*0.5);
          //vector<double> vmus(o.mus, o.mus + o.degree);
          //print(vmus);
        }


        o.L = Lfunc_init(o.degree, uint64_t(o.conductor), o.symdegree * 0.5, o.mus, &o.ecode);
        ((Lfunc *) o.L)->self_dual = YES;
        break;
      case 4:
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
  for(auto &elt: C.special_values)
    acb_clear(&elt);
  delete[] C.mus;
  Lfunc_clear(C.L);
}

// implemented at the bottom
// sets L to Sym^d L
void sympow_ECQ(fmpz_polyxx& L, const int &symdegree);

int smalljac_callback(
     __attribute__((unused)) smalljac_curve_t c,
    unsigned long q,		// prime (or prime power) in [start,end]
    int good,						// 1 if good reduction, 0 if bad
    long a[],						// n coefficients of L_q(T) with a[0] = a_1
    int n,							// either 1 or g for good p, 0 for bad p
    void *arg){  				// forwarded arg from caller

	curve *C = (curve *)arg;
  static acb_poly_t local_factor;
  static fmpz_polyxx local_factor_zz;
  static bool init = false;
  if(!init) {
    acb_poly_init(local_factor);
    local_factor_zz.fit_length(MAX_DEGREE + 1);
    acb_poly_fit_length(local_factor, MAX_DEGREE + 1);
    init = true;
  }
  local_factor_zz = 0;
  if( good ) {
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
    for(size_t i = 0; i < it->second.size(); ++i)
      local_factor_zz.set_coeff(i, it->second[i]);
    C->bad_factors.erase(it); // we will no longer use it
  }
  bool use_lpoly = true;
  if(C->nf_degree == 2){
    if(n_is_square(q)) {
        q = n_sqrt(q); // so we call Lfunc_use_lpoly accordingly
        fmpz_polyxx local_factor_zz2(2*local_factor_zz.degree() + 1);
        for(long i=0; i <= 2*local_factor_zz.degree(); ++i) {
          if(i%2 == 0)
            local_factor_zz2.set_coeff(i, local_factor_zz.get_coeff(i/2));
        }
        local_factor_zz = local_factor_zz2;
    } else { // p is split, as smalljac doesn't handle ramified primes in the monic order
      if(C->lastq == q) {
        local_factor_zz *= C->last_local_factor;
      } else {
        C->last_local_factor = fmpz_polyxx(local_factor_zz);
        use_lpoly = false;
      }
    }
    C->lastq = q;
  }
  if(use_lpoly) {
    if(C->symdegree > 1 and good)
      sympow_ECQ(local_factor_zz, C->symdegree);

    _acb_poly_set_length(local_factor, local_factor_zz.degree() + 1);
    for(long i = 0; i <= local_factor_zz.degree(); ++i) {
      // acb_poly_get_coeff_ptr(local_factor, i) = local_factor->coeffs + i
      acb_set_fmpz(local_factor->coeffs + i, local_factor_zz.coeff(i)._fmpz());
    }
    Lfunc_use_lpoly(C->L, q, local_factor);
    C->touched[q] = true;
    if(local_factor_zz.degree() >= 1) {
      if( q <= APHASH_MAXP ) {
        C->ap[primes_dict[(int16_t)q]] = -local_factor_zz.coeff(1);
      }
      // += ap^2 / p^w = (ap/n^(w/2))^2
      C->second_moment += pow(local_factor_zz.coeff(1).to<double>(), 2)/pow(double(q), C->symdegree);
    }
  }
  //acb_poly_clear(local_factor);
  return true;
}

long ap(long p, void *arg) {
	curve *C = (curve *)arg;
  return (C->ap[primes_dict[p]] % M61UI).to<slong>();
}

// what does this return?
long populate_local_factors(curve &C) {
  size_t bound = Lfunc_nmax(C.L) < 1<<13 ? 1<<13 : Lfunc_nmax(C.L);
  C.touched = vector<bool>(bound + 1, false);
  long res = smalljac_Lpolys(C.sj_curve, 1, bound, 0, smalljac_callback, &C);
  if(C.nf_degree == 2) { // deal with ramified primes
    acb_poly_t local_factor;
    acb_poly_init(local_factor);
    acb_poly_fit_length(local_factor, C.degree + 1);
    for(const auto &it : C.bad_factors){
      assert_print(C.nf_disc % fmpzxx(it.first), ==, 0);
      if( it.first <= (long) bound ) {
        _acb_poly_set_length(local_factor, it.second.size());
        for(size_t i = 0; i < it.second.size(); ++i) {
          // acb_poly_get_coeff_ptr(local_factor, i) = local_factor->coeffs + i
          acb_set_fmpz(local_factor->coeffs + i, it.second[i]._fmpz());
        }
        Lfunc_use_lpoly(C.L, it.first, local_factor);
        C.touched[it.first] = true;
        if(it.second.size() > 1) {
          if( it.first <= APHASH_MAXP )
            C.ap[primes_dict[(int16_t)it.first]] = it.second[1];
          // += ap^2 / p^w = (ap/n^(w/2))^2
          C.second_moment += pow(it.second[1].to<double>(), 2u)/ pow(double(it.first), C.symdegree);
        }
      }
    }
    // deal with inert primes
    // local_factor = 1 + ?T^2, and p^2 > bound;
    // thus we are only setting coefficients to zero
    // no need to worry about second_moment or ap
    acb_poly_one(local_factor);
    primesieve_iterator it;
    primesieve_init(&it);
    uint64_t p=0;
    while((p=primesieve_next_prime(&it)) <= bound){
      if( not C.touched[p] ){
        Lfunc_use_lpoly(C.L, p, local_factor);
        C.touched[p] = true;
      }
    }
    primesieve_free_iterator(&it);
    acb_poly_clear(local_factor);
  }
  C.tracehash = aphash(ap, &C);
  long primepi = 0;
  for(auto elt: C.touched) if(elt) ++primepi;
  C.second_moment /= primepi;
  return res;
}



ostream& operator<<(ostream &s, curve &C) {
  Lfunc_t &L = C.L;

  s << C.label <<":";
  // trace hash
  s << C.tracehash << ":";
  // second moment
  s << std::setprecision(17) << C.second_moment << ":";
  // root number
  s << Lfunc_epsilon(L) <<":";
  // r = rank
  s << Lfunc_rank(L) << ":";
  // L(1/2)^r / r! as arb
  s << Lfunc_Taylor(L) << ":";
   // special values
  s << C.special_values << ":";
  // first zeros as balls, the rest as doubles if we have enough precision
  ostream_zeros(s, L, 0);
  s << ":";
  Lplot_t *Lpp = Lfunc_plot_data(L, 0, 64.0/C.degree, 257);
  s << Lpp;
  Lfunc_clear_plot(Lpp);
  return s;
}




int main (int argc, char**argv) {
  vector<size_t> primes_vector(APHASH_MAXPI, 0);
  primes = (int16_t *)primesieve_generate_n_primes(APHASH_MAXPI, 0, INT16_PRIMES);
  for(size_t i=0; i < APHASH_MAXPI; ++i)
    primes_dict[primes[i]] = i;
  assert_print(primes[APHASH_MAXPI-1], ==, 8191);

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

      // we use printn to match SAGE's _repr_
      printf("\tRank = %" PRIu64 "\n",Lfunc_rank(L));
      printf("\tEpsilon = ");acb_printn(Lfunc_epsilon(L) ,20, 0);printf("\n");
      printf("\tFirst non-zero Taylor coeff = ");arb_printn(Lfunc_Taylor(L), 20, 0);printf("\n");
      printf("\tFirst zero = ");arb_printn(Lfunc_zeros(L, 0), 20, 0);printf("\n");


      double shift = C.symdegree*0.5;
      for(size_t i = 0; i < C.special_values.size(); ++i) {
        acb_init(&C.special_values[i]);
        double val = 1 + i + shift;
        ecode |= Lfunc_special_value(&C.special_values[i], L, val, 0);
        if(fatal_error(ecode)) {
          fprint_errors(stderr,ecode);
          std::abort();
        }
        printf("\tL(%.2f) = ", val);acb_printn(&C.special_values[i],20, 0);printf("\n");
      }
      //printf("\tFirst 20 zeros\n");
      //arb_srcptr zeros=Lfunc_zeros(L, 0);
      //for(int i  = 0; i < 20; ++i) {
      //  printf("\tZero %d = ", i);
      //  arb_printn(zeros+i, 20, 0);
      //  printf("\n");
      //}


      output << C << endl;
      // print any warnings collected along the way
      // ignore could not achieve target error bound in special value
      if( ecode != ERR_SUCCESS and ecode != ERR_SPEC_PREC ) {
        cerr << "\tBegin warnings for " << C.label << endl;
        fprint_errors(stderr, ecode);
        cerr << "\tEnd warnings for " << C.label << endl;
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
  primesieve_free(primes);
}



void sympow_ECQ(fmpz_polyxx& L, const int &symdegree) {
  assert_print(symdegree, <=, 8);
  assert_print(L.degree(), ==, 2);
  assert_print(L.get_coeff(0).is_one(), ==, true);
  static array<fmpzxx, 21> a;
  static array<fmpzxx, 37> p;
  const array<int, 7> maxapower = {2, 4, 6, 9, 12, 16, 20};
  const array<int, 7> maxppower = {3, 6, 10, 15, 21, 28, 36};

  static bool init = false;
  if(!init) {
    init = true;
    for(auto &elt: a)
      elt.set_zero();
    for(auto &elt: p)
      elt.set_zero();
    a[0].set_one();
    p[0].set_one();
  }
  a[1] = -L.get_coeff(1); // a
  p[1] = L.get_coeff(2); // p

  if(symdegree%2 == 1) {
    for(int i=2; i <= maxapower[symdegree-2]; ++i )
          a[i] = a[i/2]*a[i - i/2];
  } else {
    a[2] = a[1] * a[1];
    // we only need the even powers
    for(int i=4; i <= maxapower[symdegree-2]; i+=2 )
          a[i] = a[2*(i/4)]*a[i - 2*(i/4)];
  }



  // we can skip some of the larger powers
  int symdegreechalf = ceil(symdegree*0.5);
  for(int i=2; i <= maxppower[symdegree-2] - symdegreechalf; ++i ) {
    p[i] = p[i/2]*p[i - i/2];
  }
  int i = maxppower[symdegree-2];
  p[i] = p[i/2]*p[i - i/2];
  switch(symdegree) {
    case 2:
      L.fit_length(4);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[2] + p[1]);
      L.set_coeff(2, a[2]*p[1] - p[2]);
      L.set_coeff(3, -p[3]);
      break;
    case 3:
      L.fit_length(5);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[3] + 2*a[1]*p[1]);
      L.set_coeff(2, a[4]*p[1] - 3*a[2]*p[2] + 2*p[3]);
      L.set_coeff(3, -a[3]*p[3] + 2*a[1]*p[4]);
      L.set_coeff(4, p[6]);
      break;
    case 4:
      L.fit_length(6);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[4] + 3*a[2]*p[1] - p[2]);
      L.set_coeff(2, a[6]*p[1] - 5*a[4]*p[2] + 7*a[2]*p[3] - 2*p[4]);
      L.set_coeff(3, -a[6]*p[3] + 5*a[4]*p[4] - 7*a[2]*p[5] + 2*p[6]);
      L.set_coeff(4, a[4]*p[6] - 3*a[2]*p[7] + p[8]);
      L.set_coeff(5, -p[10]);
      break;
    case 5:
      L.fit_length(7);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[5] + 4*a[3]*p[1] - 3*a[1]*p[2]);
      L.set_coeff(2, a[8]*p[1] - 7*a[6]*p[2] + 16*a[4]*p[3] - 13*a[2]*p[4] + 3*p[5]);
      L.set_coeff(3, -a[9]*p[3] + 8*a[7]*p[4] - 22*a[5]*p[5] + 23*a[3]*p[6] - 6*a[1]*p[7]);
      L.set_coeff(4, a[8]*p[6] - 7*a[6]*p[7] + 16*a[4]*p[8] - 13*a[2]*p[9] + 3*p[10]);
      L.set_coeff(5, -a[5]*p[10] + 4*a[3]*p[11] - 3*a[1]*p[12]);
      L.set_coeff(6, p[15]);
      break;
    case 6:
      L.fit_length(8);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[6] + 5*a[4]*p[1] - 6*a[2]*p[2] + p[3]);
      L.set_coeff(2, a[10]*p[1] - 9*a[8]*p[2] + 29*a[6]*p[3] - 40*a[4]*p[4] + 22*a[2]*p[5] - 3*p[6]);
      L.set_coeff(3, -a[12]*p[3] + 11*a[10]*p[4] - 46*a[8]*p[5] + 90*a[6]*p[6] - 81*a[4]*p[7] + 28*a[2]*p[8] - 3*p[9]);
      L.set_coeff(4, a[12]*p[6] - 11*a[10]*p[7] + 46*a[8]*p[8] - 90*a[6]*p[9] + 81*a[4]*p[10] - 28*a[2]*p[11] + 3*p[12]);
      L.set_coeff(5, -a[10]*p[10] + 9*a[8]*p[11] - 29*a[6]*p[12] + 40*a[4]*p[13] - 22*a[2]*p[14] + 3*p[15]);
      L.set_coeff(6, a[6]*p[15] - 5*a[4]*p[16] + 6*a[2]*p[17] - p[18]);
      L.set_coeff(7, -p[21]);
      break;
    case 7:
      L.fit_length(9);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[7] + 6*a[5]*p[1] - 10*a[3]*p[2] + 4*a[1]*p[3]);
      L.set_coeff(2, a[12]*p[1] - 11*a[10]*p[2] + 46*a[8]*p[3] - 91*a[6]*p[4] + 86*a[4]*p[5] - 34*a[2]*p[6] + 4*p[7]);
      L.set_coeff(3, -a[15]*p[3] + 14*a[13]*p[4] - 79*a[11]*p[5] + 229*a[9]*p[6] - 359*a[7]*p[7] + 292*a[5]*p[8] - 106*a[3]*p[9] + 12*a[1]*p[10]);
      L.set_coeff(4, a[16]*p[6] - 15*a[14]*p[7] + 92*a[12]*p[8] - 296*a[10]*p[9] + 533*a[8]*p[10] - 532*a[6]*p[11] + 277*a[4]*p[12] - 68*a[2]*p[13] + 6*p[14]);
      L.set_coeff(5, -a[15]*p[10] + 14*a[13]*p[11] - 79*a[11]*p[12] + 229*a[9]*p[13] - 359*a[7]*p[14] + 292*a[5]*p[15] - 106*a[3]*p[16] + 12*a[1]*p[17]);
      L.set_coeff(6, a[12]*p[15] - 11*a[10]*p[16] + 46*a[8]*p[17] - 91*a[6]*p[18] + 86*a[4]*p[19] - 34*a[2]*p[20] + 4*p[21]);
      L.set_coeff(7, -a[7]*p[21] + 6*a[5]*p[22] - 10*a[3]*p[23] + 4*a[1]*p[24]);
      L.set_coeff(8, p[28]);
      break;
    case 8:
      L.fit_length(10);
      L.set_coeff(0, 1);
      L.set_coeff(1, -a[8] + 7*a[6]*p[1] - 15*a[4]*p[2] + 10*a[2]*p[3] - p[4]);
      L.set_coeff(2, a[14]*p[1] - 13*a[12]*p[2] + 67*a[10]*p[3] - 174*a[8]*p[4] + 239*a[6]*p[5] - 166*a[4]*p[6] + 50*a[2]*p[7] - 4*p[8]);
      L.set_coeff(3, -a[18]*p[3] + 17*a[16]*p[4] - 121*a[14]*p[5] + 467*a[12]*p[6] - 1057*a[10]*p[7] + 1415*a[8]*p[8] - 1073*a[6]*p[9] + 416*a[4]*p[10] - 70*a[2]*p[11] + 4*p[12]);
      L.set_coeff(4, a[20]*p[6] - 19*a[18]*p[7] + 154*a[16]*p[8] - 694*a[14]*p[9] + 1900*a[12]*p[10] - 3244*a[10]*p[11] + 3416*a[8]*p[12] - 2121*a[6]*p[13] + 711*a[4]*p[14] - 110*a[2]*p[15] + 6*p[16]);
      L.set_coeff(5, -a[20]*p[10] + 19*a[18]*p[11] - 154*a[16]*p[12] + 694*a[14]*p[13] - 1900*a[12]*p[14] + 3244*a[10]*p[15] - 3416*a[8]*p[16] + 2121*a[6]*p[17] - 711*a[4]*p[18] + 110*a[2]*p[19] - 6*p[20]);
      L.set_coeff(6, a[18]*p[15] - 17*a[16]*p[16] + 121*a[14]*p[17] - 467*a[12]*p[18] + 1057*a[10]*p[19] - 1415*a[8]*p[20] + 1073*a[6]*p[21] - 416*a[4]*p[22] + 70*a[2]*p[23] - 4*p[24]);
      L.set_coeff(7, -a[14]*p[21] + 13*a[12]*p[22] - 67*a[10]*p[23] + 174*a[8]*p[24] - 239*a[6]*p[25] + 166*a[4]*p[26] - 50*a[2]*p[27] + 4*p[28]);
      L.set_coeff(8, a[8]*p[28] - 7*a[6]*p[29] + 15*a[4]*p[30] - 10*a[2]*p[31] + p[32]);
      L.set_coeff(9, -p[36]);
      break;
    default:
      throw_line("we cannot get here!"s);
  }
}


/* You may generate the function above in sage with:
def tensor_charpoly(f, g):
    R = PolynomialRing(g.parent(), "y");
    y = R.gen();
    A = f(y)
    B = R(g.homogenize(y))
    return B.resultant(A)
def base_change(Lpoly, r):
    R = Lpoly.parent()
    T = R.gen()
    S = PolynomialRing(R, 'u')
    u = S.gen()
    return R(Lpoly(u).resultant(u**r - T))
def sym_pol(L, n):
    Ln = L
    for i in range(2, n + 1):
        Ln = tensor_charpoly(Ln, L)
        b = base_change(L, i)
        extra = (Ln//b).sqrt(2)
        Ln = b*extra
        if Ln[0] == -1:
            Ln *= -1
    return Ln
S.<a, p> = ZZ[]
R.<T> = S[]
L = 1 - a*T + p*T^2


max_symdeg = 8
Ls = {i: sym_pol(L, i) for i in range(2, max_symdeg+1)}

out = ""


maxapower = []
maxppower = []
for i in range(2, max_symdeg + 1):
    maxa, maxp = [max([max([elt1.degree(var) for elt1 in elt.monomials()]) for elt in Ls[i].list()[1:]]) for var in [a, p]]
    maxapower.append(maxa)
    maxppower.append(maxp)
apowers = maxa + 1
ppowers = maxp + 1
out += r"""
assert_print(symdegree, <=, {max_symdeg});
assert_print(L.degree(), ==, 2);
assert_print(L.get_coeff(0).is_one(), ==, true);
static std::array<fmpzxx, {maxa1}> a;
static std::array<fmpzxx, {maxp1}> p;
const std::array<int, {length}> maxapower = {maxapower};
const std::array<int, {length}> maxppower = {maxppower};
""".format(max_symdeg=max_symdeg,
           maxa1=maxa + 1,
           maxp1=maxp + 1,
           maxapower = str(maxapower).replace('[','{').replace(']','}'),
           maxppower = str(maxppower).replace('[','{').replace(']','}'),
           length = len(maxapower)
          );
out += r"""
static bool init = false;
if(!init) {
  init = true;
  for(auto &elt: a)
    elt.set_zero();
  for(auto &elt: p)
    elt.set_zero();
  a[0].set_one();
  p[0].set_one();
}
a[1] = -L.get_coeff(1); // a
p[1] = L.get_coeff(2); // p

if(symdegree%2 == 1) {
  for(int i=2; i <= maxapower[symdegree-2]; ++i )
        a[i] = a[i/2]*a[i - i/2];
} else {
  a[2] = a[1] * a[1];
  // we only need the even powers
  for(int i=4; i <= maxapower[symdegree-2]; i+=2 )
        a[i] = a[2*(i/4)]*a[i - 2*(i/4)];
}

// we can skip some of the larger powers
int symdegreechalf = ceil(symdegree*0.5);
for(int i=2; i <= maxppower[symdegree-2] - symdegreechalf; ++i ) {
    p[i] = p[i/2]*p[i - i/2];
}
int i = maxppower[symdegree-2];
p[i] = p[i/2]*p[i - i/2];
switch(symdegree) {
"""

out = '\n'.join(['  '  + elt for elt in out.split('\n')])


for i in range(2, max_symdeg + 1):
    localout = ""
    localout += 'case %d:\n' % i;
    localout += '  L.fit_length(%d);\n' % (Ls[i].degree() + 1)
    for j, c in enumerate(Ls[i].list()):
        if c != 0:
            expr = str(c)
            # x -> x[1]
            expr = re.sub(r"([ap])(?!\^)", r"\1[1]", expr)
            # x^e -> x[e]
            expr = re.sub(r"([ap])\^(\d+)", r"\1[\2]", expr)
            localout +='  L.set_coeff(%d, %s);\n' % (j, expr)
    localout += '      break;\n'
    
    out += '\n'.join(['    '  + elt for elt in localout.split('\n')])
out += r"""
    default:
      throw_line("we cannot get here!"s);
    }
"""
out = 'void sympow_ECQ(fmpz_polyxx& L, const int &symdegree) {' + out + '\n}'
print(out)
*/

#endif

