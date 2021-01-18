// Copyright Edgar Costa 2018-2019
// See LICENSE file for license details.
//
/* Computes the L-function data for the LMFDB for Artin Reps given
 * the local data to deduce the local factors from the conjugacy classes
 *
 * python function to generate input provided at the bottom with an example
 */
#define __STDC_FORMAT_MACROS
#define special_values_size 2 // implies computing L(1) ... L(special_values_size)
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <flint/fmpz.h>
#include <flint/fmpzxx.h>
#include <flint/nmod_poly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/ulong_extras.h>
#include <flint/perm.h>
#include <primesieve.hpp>
#include <acb_poly.h>
#include "glfunc.h"
#include "examples_tools.h"

using std::cerr;
using std::cout;
using std::endl;
using std::exception;
using std::int8_t;
using std::int64_t;
using std::ifstream;
using std::istream;
using std::iswspace;
using std::getline;
using std::lexicographical_compare;
using std::make_pair;
using std::map;
using std::ostream;
using std::ofstream;
using std::pair;
using std::set;
using std::size_t;
using std::string;
using std::stringstream;
using std::vector;

typedef std::chrono::time_point<std::chrono::system_clock> SystemTime ;



/******************************************************************************
 * Handy tools for fq_nmod_poly_t
 *  - conv, to convert a vector<int64_t> to a Fq[x] polynomial
 *  - roots, to compute the roots of Fq[x] polynomial
 *  - roots, to compute the roots of ZZ[x] polynomial over Fq
 *****************************************************************************/



// converts vector<fmpzxx> to a fq_nmod_poly_t
// if neg, returns f(-x) or -f(-x) to keep it monic, in case it was monic
void conv(fq_nmod_poly_t res,
          const vector<fmpzxx> &polynomial,
          const fq_nmod_ctx_t ctx,
          const bool &neg = false) {
  fq_nmod_poly_fit_length(res, polynomial.size(), ctx);
  fq_nmod_t tmp;
  fq_nmod_init2(tmp, ctx);
  size_t parity = polynomial.size() % 2;

  for(size_t i = 0; i < polynomial.size(); ++i) {
    fq_nmod_set_fmpz(
        tmp,
        (neg and i % 2 == parity) ? (-polynomial[i]).evaluate()._fmpz() : polynomial[i]._fmpz(),
        ctx);
    fq_nmod_poly_set_coeff(res, i, tmp, ctx);
  }
  fq_nmod_clear(tmp, ctx);
}

// converts vector<fmpzxx> to a nmod_poly_t
void conv(nmod_poly_t pol,
          const vector<fmpzxx> &polynomial) {
  nmod_poly_fit_length(pol, polynomial.size());
  for(size_t i = 0; i < polynomial.size(); ++i)
    nmod_poly_set_coeff_ui(pol, i,
        (polynomial[i] % pol->mod.n).to<mp_limb_t>());
}

// computes the roots of split polynomial over Fq
void split_roots(vector<fq_nmod_struct> &res,
                 const vector<fmpzxx> &polynomial,
                 const fq_nmod_ctx_t ctx) {
  // if not empty, we will leak memory
  assert_print(res.size(), ==, 0);
  res.resize(polynomial.size() - 1);
  for(auto &elt : res)
    fq_nmod_init2(&elt, ctx);

  fq_nmod_poly_factor_t fac;
  fq_nmod_poly_factor_init(fac, ctx);
  fq_nmod_poly_factor_fit_length(fac, polynomial.size(), ctx);

  fq_nmod_poly_t f;
  fq_nmod_poly_init2(f, polynomial.size(), ctx);
  // f = f(-x) or -f(-x)
  conv(f, polynomial, ctx, true);

  assert(fq_nmod_is_one(f->coeffs + polynomial.size() - 1, ctx));
  fq_nmod_poly_factor_equal_deg(fac, f, 1, ctx);
  assert_print(fac->num, ==, mp_limb_signed_t(polynomial.size()) - 1);

  for(mp_limb_signed_t i = 0; i < fac->num; ++i) {
    // because we factored f(-x) or -f(-x)
    // we don't need to negate the constant coefficient
    fq_nmod_poly_get_coeff(&res[i], fac->poly + i, 0, ctx);
  }

  fq_nmod_poly_factor_clear(fac, ctx);
  fq_nmod_poly_clear(f, ctx);
}



/******************************************************************************
 * Computing Euler factors of of Artin L-functions
 *
 * This is heavily inspired by a conversion of inspired:
 * This is basically a C++ version of
 *  - lmfdb/artin_representations/cyc_alt_res_engine.py
 * originally written by Paul-Olivier Dehaye.
 *
 *
 * This provides an interface to match each prime, to a conjugacy class, and
 * hence with the corresponding Euler factor.
 *
 * There are four cases depending on p:
 *
 *  1- "hard", if p divides the polynomial of the discriminant, a finite
 *     superset of the bad primes. The list of hard primes and the
 *     corresponding Euler factors is part of the input.
 *
 *  2- "CYC", when the cycle type of Frobenius is enough to determine the conjugacy class.
 *
 * The two remaining cases are handled by using
 * technique introduced in [1].
 * In short, let f in ZZ[x] be the defining polynomial of the field with Galois
 * group G. There is a polynomial h in ZZ[x] and polynomials Gamma_C indexed by the conjugacy classes of G. Such that
 *    Frob_p \in C <=> Gamma_C( Tr Fp[x]/f(x) h(x) x^p ) = 0
 * In other, words we can compute an invariant alpha, and use a precomputed table
 * of polynomials to figure out which polynomial vanishes at alpha.
 * In practice, we try to take advantage of additional symmetries in G, see Remark 5.10 in [1], to reduce the degree of Gamma_C.
 *
 *  3- "RES", this is the generic case of the above.
 *
 *  4- "ALT", use Serre's trick to distinguish which of two conjugacy classes is
 *            the right one. A particular and simpler version of "RES".
 *            See Example 3.9 in [1].
 *
 * [1] -  Tim and Vladimir Dokchitser, "Identifying Frobenius elements in Galois groups".
 */


// given p, figures out the cycle type of Frobenious
void cycle_type(vector<size_t>& res,
                const vector<fmpzxx> &polynomial,
                const int64_t &p) {
  nmod_poly_t pol;
  nmod_poly_init2(pol, p, polynomial.size());
  conv(pol, polynomial);
  nmod_poly_factor_t fac;
  nmod_poly_factor_init(fac);
  nmod_poly_factor(fac, pol);

  res = vector<size_t>(fac->num, 0);
  for(mp_limb_signed_t i = 0; i < fac->num; ++i){
    res[i] = nmod_poly_degree(fac->p + i);
    assert_print(fac->exp[i], ==, 1);
  }
  sort(res.begin(), res.end());
  nmod_poly_factor_clear(fac);
  nmod_poly_clear(pol);
}

//returns sum([r ** (p ** j) for j in powers])
void powers_sum(fq_nmod_t res,
                const fq_nmod_t r,
                const vector<int64_t> &powers,
                const fq_nmod_ctx_t ctx) {
  fq_nmod_zero(res, ctx);
  const fmpz *p = fq_nmod_ctx_prime(ctx);
  fmpz_t tmp;
  fmpz_init(tmp);
  fq_nmod_t power;
  fq_nmod_init2(power, ctx);
  for(auto &j : powers) {
    fmpz_pow_ui(tmp, p, j);
    fq_nmod_pow(power, r, tmp, ctx);
    fq_nmod_add(res, res, power, ctx);
  }
  fq_nmod_clear(power, ctx);
  fmpz_clear(tmp);
}

//sum(gamma_polynomial(r) * powers_sum(r, powers, p) for r in roots)
void alpha_res(mp_limb_t &res,
               const vector<fq_nmod_struct> &roots,
               const vector<int64_t> &powers,
               const vector<fmpzxx> &gamma_polynomial,
               const fq_nmod_ctx_t ctx) {

  fq_nmod_poly_t gamma;
  fq_nmod_poly_init2(gamma, gamma_polynomial.size(), ctx);
  conv(gamma, gamma_polynomial, ctx);

  fq_nmod_t sum, tmp, gr, pr;
  fq_nmod_init2(sum, ctx);
  fq_nmod_init2(gr, ctx);
  fq_nmod_init2(pr, ctx);
  fq_nmod_init2(tmp, ctx);
  fq_nmod_zero(sum, ctx);
  for(const auto &r: roots) {
    powers_sum(pr, &r, powers, ctx);
    fq_nmod_poly_evaluate_fq_nmod(gr, gamma, &r, ctx);
    fq_nmod_mul(tmp, pr, gr, ctx);
    fq_nmod_add(sum, sum, tmp, ctx);
  }
  fq_nmod_clear(gr, ctx);
  fq_nmod_clear(pr, ctx);
  fq_nmod_clear(tmp, ctx);
  assert_print(nmod_poly_length(sum), <=, 1);
  if(nmod_poly_length(sum) == 0)
    res = 0;
  else
    res = sum->coeffs[0];
  fq_nmod_clear(sum, ctx);
  fq_nmod_poly_clear(gamma, ctx);
}





void conjugacy_class_matcher_res_alt(size_t &c,
                                     const mp_limb_t &alpha,
                                     const vector< pair<size_t, vector<fmpzxx> > > *root_of,
                                     const int64_t &p) {
  size_t max_len = 0;
  for(const auto &elt: *root_of) {
    if(elt.second.size() > max_len)
      max_len = elt.second.size();
  }

  c = SIZE_MAX;

  nmod_poly_t gamma_c;
  nmod_poly_init2(gamma_c, p, max_len);
  for(const auto &elt: *root_of) {
    conv(gamma_c, elt.second);
    if( 0 == nmod_poly_evaluate_nmod(gamma_c, alpha) ) {
      c = elt.first;
      break;
    }
  }
  assert_print(c, !=, SIZE_MAX);
  nmod_poly_clear(gamma_c);
}


/*
 * Some basic code to figure out if two permutation are conjugate in A_n
 */

// assumes perm[i] \in [0..perm.size() - 1]
void cycle_tuples(vector< vector<size_t> > &res, const vector<size_t> &perm) {
  const size_t &n = perm.size();
  res.resize(0);
  vector<bool> seen(n, false);
  for(size_t i = 0; i < n; ++i) {
    if( not seen[i] ) {
      if(perm[i] == i) {
        res.push_back(vector<size_t>(1, i));
        seen[i] = true;
      } else {
        vector<size_t> cycle(1, i);
        size_t k = perm[i];
        while( k != i ) {
          cycle.push_back(k);
          seen[k] = true;
          k = perm[k];
        }
        res.push_back(cycle);
      }
    }
  }
  // sort using the length of the cycles
  sort(res.begin(),
      res.end(),
      [](const vector<size_t> &a, const vector<size_t> &b) {
      return a.size() > b.size();
      });
}

void cycle_types(vector<size_t>& res, const vector< vector<size_t>> cycle_tuples) {
  res = vector<size_t>(cycle_tuples.size(), 0);
  for(size_t i = 0; i < res.size(); ++i)
    res[i] = cycle_tuples[i].size();
}


// assumes perm[i] \in [0..perm.size() - 1] for perm = rho or sigma
bool are_conjugate_An(const vector<size_t> &rho, const vector<size_t> &sigma) {
  vector< vector<size_t> > rho_cycle_tuple, sigma_cycle_tuple;
  cycle_tuples(rho_cycle_tuple, rho);
  cycle_tuples(sigma_cycle_tuple, sigma);

  vector<size_t> rho_cycle_type, sigma_cycle_type;
  cycle_types(rho_cycle_type, rho_cycle_tuple);
  cycle_types(sigma_cycle_type, sigma_cycle_tuple);

  if ( rho_cycle_type !=  sigma_cycle_type )
    return false;

  // Conjugacy classes in A_n are the same as in S_n,
  // except when all the cycles are of odd and different lengths.
  // In that case, S_n conjugacy classes split up in 2 A_n conjugacy classes
  if ( set<size_t>(sigma_cycle_type.begin(), sigma_cycle_type.end()).size()
      != sigma_cycle_type.size() )
    return true;
  size_t n = 0;
  for(auto &x : sigma_cycle_type) {
    if ( x % 2 == 0 )
      return true;
    n += x;
  }

  // Now we know we only have odd length cycles, all different
  // cycles match up with their lengths, as the tuples were already sorted
  slong * tmp = new slong[n];
  for(size_t i = 0; i < n; ++i)
    tmp[i] = i;
  for(size_t i = 0; i < rho_cycle_tuple.size(); ++i) {
      const vector<size_t> &rho_tuple = rho_cycle_tuple[i];
      const vector<size_t> &sigma_tuple = sigma_cycle_tuple[i];
      for(size_t j = 0; j < rho_tuple.size(); ++j)
          tmp[rho_tuple[j]] = sigma_tuple[j];
  }
  int sig = _perm_parity(tmp, n);
  delete[] tmp;
  return sig == 0;
}


inline bool operator<(const pair<size_t, const fq_nmod_struct*>& a,
                      const pair<size_t, const fq_nmod_struct*>& b) {
  return lexicographical_compare(
      a.second->coeffs, a.second->coeffs + a.second->length,
      b.second->coeffs, b.second->coeffs + b.second->length);
}


void frobenius_permutation(vector<size_t> &res,
                           const vector<fq_nmod_struct> &roots,
                           const fq_nmod_ctx_t ctx) {
  uint64_t p = fmpz_get_ui(fq_nmod_ctx_prime(ctx));
  vector<pair<size_t, const fq_nmod_struct*> > roots_ext(roots.size());
  vector<pair<size_t, const fq_nmod_struct*> >roots_frob_ext(roots.size());
  vector<fq_nmod_struct> roots_frob(roots.size());
  for(size_t i = 0; i < roots.size(); ++i) {
    roots_ext[i] = make_pair(i, &roots[i]);

    fq_nmod_init2(&roots_frob[i], ctx);
    fq_nmod_pow_ui(&roots_frob[i], &roots[i], p, ctx);
    roots_frob_ext[i] = make_pair(i, &roots_frob[i]);
  }

  sort(roots_ext.begin(), roots_ext.end());
  sort(roots_frob_ext.begin(), roots_frob_ext.end());
  res = vector<size_t>(roots.size());
  for(size_t i = 0; i < roots.size(); ++i) {
    res[roots_frob_ext[i].first] = roots_ext[i].first;
  }
  for(auto &elt: roots_frob)
    fq_nmod_clear(&elt, ctx);
}



void alpha_alt(mp_limb_t &res,
               const vector<fq_nmod_struct> &roots,
               const vector<size_t> &permutation,
               const fq_nmod_ctx_t ctx) {
  // see Example 3.9 in [1]
  vector<size_t> frob_perm;
  frobenius_permutation(frob_perm, roots, ctx);
  // prod = prod_{i < j} (roots[i] - roots[j])
  fq_nmod_t prod, diff;
  fq_nmod_init2(diff, ctx);
  fq_nmod_init2(prod, ctx);
  fq_nmod_one(prod, ctx);
  for(size_t i = 0; i < roots.size(); ++i) {
    for(size_t j = 0; j < i; ++j) {
      fq_nmod_sub(diff, &roots[i], &roots[j], ctx);
      fq_nmod_mul(prod, prod, diff, ctx);
    }
  }
  assert_print(nmod_poly_length(prod), ==, 1);
  res = prod->coeffs[0];
  fq_nmod_clear(diff, ctx);
  fq_nmod_clear(prod, ctx);
  if( not are_conjugate_An(frob_perm, permutation) )
    res = nmod_neg(res, ctx->mod);
}


void conv(vector<acb_poly_struct> &local_factors, const vector< vector< vector< fmpzxx > > > &local_factors_int, const size_t &n, const int64_t &prec) {
  assert_print(local_factors.size(), ==, 0);
  local_factors.resize(local_factors_int.size());
  acb_t tmp;
  acb_init(tmp);

  vector<acb_struct> z(n);
  for(auto &elt : z)
    acb_init(&elt);
  acb_one(&z[0]);
  if (n > 1) {
    // z[1] = 2/n
    acb_set_ui(&z[1], 2);
    acb_div_ui(&z[1], &z[1], n, prec);
    // z[1] = exp(2 pi i /n)
    acb_exp_pi_i(&z[1], &z[1], prec);
    for(size_t i = 2; i < n; ++i)
      acb_mul(&z[i], &z[i-1], &z[1], prec);
  }

  for(size_t i = 0; i < local_factors_int.size(); ++i) {
    const vector< vector<fmpzxx> > &pzz = local_factors_int[i];
    acb_poly_init(&local_factors[i]);
    acb_poly_fit_length(&local_factors[i], pzz.size());
    for(size_t j = 0; j < pzz.size(); ++j) {
      acb_zero(tmp);
      for(size_t k = 0; k < pzz[j].size(); ++k) {
        acb_addmul_fmpz(tmp, &z[k], (pzz[j][k])._fmpz(), prec);
      }
      acb_poly_set_coeff_acb(&local_factors[i], j, tmp);
    }
  }
  for(auto &elt : z)
    acb_clear(&elt);
  acb_clear(tmp);

}



typedef struct {
  // index(C) --> gamma_C
  vector< pair<size_t, vector<fmpzxx> > > root_of;
  vector<size_t> permutation;
} alt_data;

istream & operator>>(istream& s, alt_data& out) {
  pair<vector< pair<size_t, vector<fmpzxx> > >, vector<size_t>> p;
  if( !(s >> p) )
    throw_line("bad alt data"s);
  out.root_of = p.first;
  out.permutation = p.second;
  return s;
}

ostream& operator<<(ostream& s, const alt_data& a)
{
  s <<"{\nroot_of : ";
  s << a.root_of << endl;
  s << "permutation : ";
  s << a.permutation <<'}';
  return s;
}

typedef struct {
  // index(C) --> gamma_C
  vector< pair<size_t, vector<fmpzxx> > > root_of;
  vector<int64_t> powers;
  vector<fmpzxx> gamma_polynomial;
} res_data;

istream & operator>>(istream& s, res_data& out) {
  pair<vector< pair<size_t, vector<fmpzxx> > >, vector<vector<fmpzxx> > > p;
  if( !(s >> p) )
    throw_line("bad res data"s);
  out.root_of = p.first;
  assert_print(p.second.size(), ==, 2);
  out.powers.resize(p.second[0].size());
  for(size_t i = 0; i < p.second[0].size(); ++i) {
    out.powers[i] = p.second[0][i].to<slong>();
    assert_print(out.powers[i], >, 0);
  }
  out.gamma_polynomial = p.second[1];
  return s;
}

ostream& operator<<(ostream& s, const res_data& a)
{
  s <<"{\nroot_of : ";
  s << a.root_of << endl;
  s << "powers : ";
  s << a.powers <<endl;
  s << "gamma_polynomial : ";
  s << a.gamma_polynomial <<'}';
  return s;
}

typedef struct {
  // artin rep label
  string label;

  // dimension
  int64_t dimension;

  // conductor
  int64_t conductor;

  // artin field
  vector<fmpzxx> polynomial;

  // mus
  double* mus;

  // L-function
  Lerror_t ecode;
  Lfunc_t L;

  // the local factors are written in terms of zeta_n
  int64_t n;

  // conjugacy class -> euler factor
  vector<acb_poly_struct> local_factors;

  //  k -> L(k+1)
  array<acb_struct, special_values_size> special_values;

  // p -> conjugacy class
  map<fmpzxx, size_t> hard_primes;

  // cycle type -> type
  //  0 = cyc
  //  1 = alt
  //  2 = res
  map< vector<size_t>, size_t> type;

  // cycle type -> conjugacy class
  map< vector<size_t>, size_t> cyc;

  // cycle type -> alt/res_data
  map< vector<size_t>, alt_data> alt;
  map< vector<size_t>, res_data> res;

} artin_rep;


void artin_rep_clear(artin_rep &d) {
  for(auto &elt: d.local_factors)
    acb_poly_clear(&elt);
  for(auto &elt: d.special_values)
    acb_clear(&elt);

  delete[] d.mus;
  Lfunc_clear(d.L);
}

template <class T, class R>
void decrease_keys(map<T, R> &m) {
  map<T, R> buf;
  for(const auto& elt: m)
    buf[elt.first - 1] = elt.second;
  m = buf;
}
template <class T, class R>
void decrease_keys(vector< pair<T, R>> &m) {
  for(auto& elt: m)
    elt.first--;
}



void artin_rep_fix_indexes(artin_rep &d) {
  // the index of the inputs start at 1, but we want them to start at 0
  for(auto& elt : d.hard_primes)
    elt.second--;
  for(auto& elt: d.cyc)
    elt.second--;
  for(auto& elt: d.alt) {
    for(auto& z: elt.second.permutation)
      z--;
    decrease_keys(elt.second.root_of);
  }
  for(auto& elt: d.res)
    decrease_keys(elt.second.root_of);
}


istream & operator>>(istream & is, artin_rep &o)
{
  for(size_t i = 0; i < 11 ; ++i) {
    string s;
    if(not getline(is, s, ':'))
      throw_line("missing field"s);
    stringstream ss(s);

    //print(i);
    //print(ss.str());
    switch(i) {
      case 0:
        if(!(ss >> o.label))
          throw_line("bad label"s);
        break;
      case 1:
        if(!(ss >> o.dimension))
          throw_line("bad dimension"s);
        break;
      case 2:
        if(!(ss >> o.conductor))
          throw_line("bad conductor"s);
        break;
      case 3:
        if(!(ss >> o.polynomial))
          throw_line("bad polynomial"s);
        break;
      case 4:
        {
          vector<double> buf;
          if(!(ss >> buf))
            throw_line("bad mus"s);
          assert_print(buf.size(), ==, size_t(o.dimension));
          o.mus = new double[o.dimension];
          for(size_t i = 0; i < buf.size(); ++i)
            o.mus[i] = buf[i];
          // alg = anal so normalisation = 0.0
          o.L = Lfunc_init(o.dimension, o.conductor, 0.0, o.mus, &o.ecode);
          if(fatal_error(o.ecode)) {
            fprint_errors(stderr, o.ecode);
            throw_line("could not init L-function"s);
          }
        }
        break;
      case 5:
        if(!(ss >> o.n))
          throw_line("bad n"s);
        break;
      case 6:
        {
          vector< vector< vector<fmpzxx> > > local_factors_int;
          if(!(ss >> local_factors_int))
            throw_line("bad local factors"s);
          conv(o.local_factors, local_factors_int, o.n, Lfunc_wprec(o.L));
          break;
        }
      case 7:
        if(!(ss >> o.hard_primes))
          throw_line("bad hard primes"s);
        break;
      case 8:
        if(!(ss >> o.cyc))
          throw_line("bad cyc"s);
        break;
      case 9:
        if(!(ss >> o.alt))
          throw_line("bad alt"s);
        break;
      case 10:
        if(!(ss >> o.res))
          throw_line("bad res"s);
        break;
      default:
        throw_line("too many fields in the line!"s);
    }
  }
  artin_rep_fix_indexes(o);
  // populate types
  for(const auto& elt: o.cyc)
    o.type[elt.first] = 0;
  for(const auto& elt: o.alt)
    o.type[elt.first] = 1;
  for(const auto& elt: o.res)
    o.type[elt.first] = 2;
  return is;
}




ostream& operator<<(ostream &s, artin_rep &AR) {
  Lfunc_t &L = AR.L;

  s << AR.label <<":";
  // root number
  s << Lfunc_epsilon(L) <<":";
  // r = rank
  s << Lfunc_rank(L) << ":";
  // L(1/2)^r / r! as arb
  s << Lfunc_Taylor(L) << ":";
  // special values
  s << AR.special_values << ":";
  // first zeros as balls, the rest as doubles if we have enough precision
  ostream_zeros(s, L, 0);
  s << ":";
  Lplot_t *Lpp = Lfunc_plot_data(L, 0, 64.0/AR.dimension, 257);
  s << Lpp;
  Lfunc_clear_plot(Lpp);
  return s;
}



// compute the Euler poly for p
void lpoly(const int64_t &p, artin_rep &AR) {
  size_t c;
  // is a hard prime?
  auto hp = AR.hard_primes.find(fmpzxx(p));
  if( hp != AR.hard_primes.end() ) {
    c = hp->second;
  } else {
    // compute cycle type
    vector<size_t> ct;
    cycle_type(ct, AR.polynomial, p);
    size_t type = AR.type[ct];
    if( type == 0){
      c = AR.cyc[ct];
    } else {
      slong d = slong(lcm(ct));
      fq_nmod_ctx_t ctx;
      fmpz_t pz;
      fmpz_init(pz);
      fmpz_set_ui(pz, p);
      fq_nmod_ctx_init(ctx, pz, d, "a");
      fmpz_clear(pz);

      vector<fq_nmod_struct> roots;
      split_roots(roots, AR.polynomial, ctx);

      mp_limb_t alpha;
      const vector< pair<size_t, vector<fmpzxx> > > *root_of_data;
      if(type == 2) {
        const res_data &data = AR.res[ct];
        root_of_data = &data.root_of;
        alpha_res(alpha, roots, data.powers, data.gamma_polynomial, ctx);
      } else {
        const alt_data &data = AR.alt[ct];
        root_of_data = &data.root_of;
        alpha_alt(alpha, roots, data.permutation, ctx);
      }
      conjugacy_class_matcher_res_alt(c, alpha, root_of_data, p);
      for(auto &elt : roots)
        fq_nmod_clear(&elt, ctx);
      fq_nmod_ctx_clear(ctx);
    }
  }
  Lfunc_use_lpoly(AR.L, p, &AR.local_factors[c]);
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

    primesieve::iterator ps;

    int r = 0;


    while(std::getline(input, line)) {
      SystemTime start(std::chrono::system_clock::now());
      std::time_t startt = std::chrono::system_clock::to_time_t(start);
      cout << "Date:   \t" <<  std::put_time(std::localtime(&startt), "%F %T") << endl;

      artin_rep AR;
      Lerror_t &ecode = AR.ecode;
      Lfunc_t &L = AR.L;


      // read a line
      stringstream linestream(line);
      linestream >> AR;
      cout << "Starting:\t"<<AR.label<<endl;

      // we need all the local factors p <= target_M
      uint64_t target_M = Lfunc_nmax(L);

      cout <<"using p <= " << target_M << endl;

      // populate local factors
      ps.skipto(0);
      for(uint64_t p = ps.next_prime(); p <= target_M; p = ps.next_prime())
        lpoly(int64_t(p), AR);


      // do the computation
      ecode |= Lfunc_compute(L);
      if(fatal_error(ecode)) {
        fprint_errors(stderr, ecode);
        std::abort();
      }

      printf("Rank = %" PRIu64 "\n",Lfunc_rank(L));
      printf("Epsilon = ");acb_printd(Lfunc_epsilon(L),20);printf("\n");
      printf("Leading Taylor coeff = ");arb_printd(Lfunc_Taylor(L), 20);printf("\n");
      printf("First zero = ");arb_printd(Lfunc_zeros(L, 0), 20);printf("\n");

      for(size_t i = 0; i < AR.special_values.size(); ++i) {
        acb_init(&AR.special_values[i]);
        ecode |= Lfunc_special_value(&AR.special_values[i], L, i + 1, 0);
        if(fatal_error(ecode)) {
          fprint_errors(stderr,ecode);
          std::abort();
        }
        printf("L(%lu) = ", i + 1);acb_printd(&AR.special_values[i],20);printf("\n");
      }

      output << AR << endl;
      // print any warnings collected along the way
      // ignore could not achieve target error bound in special value
      if( ecode != ERR_SUCCESS and ecode != ERR_SPEC_PREC ) {
        cerr << "Begin warnings for " << AR.label << endl;
        fprint_errors(stderr, ecode);
        cerr << "End warnings for " << AR.label << endl;
        r++;
      }

      SystemTime end(std::chrono::system_clock::now());
      std::time_t endt = std::chrono::system_clock::to_time_t(end);
      double walltime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
      cout << "Date:   \t" <<  std::put_time(std::localtime(&endt), "%F %T") << endl;
      cout << "Done:   \t"<< AR.label << "\ttook ";
      cout << std::setw(6) << std::setfill(' ')  << std::fixed << std::setprecision(2) << walltime/1000 << "s"<< endl << endl;

      //free memory
      artin_rep_clear(AR);
    }
    flint_cleanup();
    return r;
  } catch( const std::exception & ex ) {
     cerr << "Uncaught exception: " <<ex.what() << endl;
     std::abort();
  }
}



/* Sample way to generate the input

from lmfdb import db
from lmfdb.utils.utilities import num2letters
import json
def json_hack(x):
    if isinstance(x,str):
        return x
    if isinstance(x, dict):
        keys = []
        values = []
        for k, v in x.items():
            keys.append(k)
            values.append(v)
        x = [keys, values]
        return json.dumps(x).replace('(','[').replace(')',']')
    return json.dumps(x)
def json_frob_resolvents(data):
    def classes(root_of):
        return [[elt['ConjugacyClass'],
                      list(map(int, elt['RootOf']))]
                        for elt in root_of]
    cyc = {}
    res = {}
    alt = {}
    for elt in data:
        alg = elt['Algorithm']
        ct = tuple(elt['CycleType'])
        if alg == 'CYC':
            cyc[ct] = elt['Classes']
        elif alg == 'ALT':
            alt[ct] = [classes(elt['Classes']), elt['Data']]
        else:
            res[ct] = [classes(elt['Classes']), [elt['Data']['Powers'],
                                                 elt['Data']['Resolvent']]]
    return ":".join([json_hack(elt) for elt in [cyc, alt, res]])
def generate_lines(AR):
    output = []
    base_label = str(AR['Baselabel'])
    dim = int(AR['Dim'])
    cond = int(AR['Conductor'])
    polynomial = list(map(int, AR['NFGal']))
    local_data = db.artin_field_data.lucky({'Polynomial':polynomial})
    n = int(AR['CharacterField'])
    hard_primes = list(map(int, AR['HardPrimes']))
    output = ""
    for elt in AR['GaloisConjugates']:
        label = base_label + "." + num2letters(elt['GalOrbIndex'])
        # is stored in DB as an index starting at
        complex_conjugation = elt['Character'][local_data['ComplexConjugation'] - 1]
        assert(complex_conjugation)
        trace_cc = complex_conjugation[0]
        mus = [0r for _ in range(int(dim + trace_cc) / 2)] + [1r for _ in range(int(dim - trace_cc) / 2)]
        local_factors = [[[int(x) for x in y] for y in z] for z in elt['LocalFactors']]
        hard_factors = dict(zip(hard_primes, map(int, elt['HardFactors'])))
        line = [label, dim, cond, polynomial, mus, n, local_factors, hard_factors]
        line_str = ":".join(map(json_hack, line))
        line_str += ":" + json_frob_resolvents(local_data['FrobResolvents'])
        output += line_str + '\n'
    return output

def export_search(query, filename):
    with open(filename, 'w') as F:
        for AR in db.artin_reps.search(query):
            F.write(generate_lines(AR))

# Example:
sage: print(generate_lines(db.artin_reps.lucky({'Conductor':47, 'Dim': 2})))
2.47.5t2.a.a:2:47:[1, 0, -1, 2, -2, 1]:[0, 1]:5:[[[1], [-2], [1]], [[1], [0], [-1]], [[1], [0, 0, -1, -1], [1]], [[1], [1, 0, 1, 1], [1]], [[1], [-1]]]:[[2, 47], [3, 5]]:[[[1, 1, 1, 1, 1], [1, 2, 2]], [1, 2]]:[[[5]], [[[[4, [-47, 1]], [3, [47, 1]]], [2, 3, 4, 5, 1]]]]:[[], []]
2.47.5t2.a.b:2:47:[1, 0, -1, 2, -2, 1]:[0, 1]:5:[[[1], [-2], [1]], [[1], [0], [-1]], [[1], [1, 0, 1, 1], [1]], [[1], [0, 0, -1, -1], [1]], [[1], [-1]]]:[[2, 47], [3, 5]]:[[[1, 1, 1, 1, 1], [1, 2, 2]], [1, 2]]:[[[5]], [[[[4, [-47, 1]], [3, [47, 1]]], [2, 3, 4, 5, 1]]]]:[[], []]

sage: print(generate_lines(db.artin_reps.lucky({'Baselabel':'3.3003175.6t11.b'})))
3.3003175.6t11.b.a:3:3003175:[1, -1, 3, 3, 3, -1, 1]:[0, 0, 1]:1:[[[1], [-3], [3], [-1]], [[1], [3], [3], [1]], [[1], [-1], [-1], [1]], [[1], [1], [-1], [-1]], [[1], [1], [-1], [-1]], [[1], [-1], [-1], [1]], [[1], [0], [0], [-1]], [[1], [1], [1], [1]], [[1], [-1], [1], [-1]], [[1], [0], [0], [1]], [[1], [-1]], [[1], [0], [-1]]]:[[3, 5, 7, 131], [8, 11, 12, 11]]:[[[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 2], [1, 1, 4], [2, 4], [3, 3], [6]], [1, 3, 8, 9, 7, 10]]:[[], []]:[[[1, 1, 2, 2], [2, 2, 2]], [[[[4, [-353, 22, 15, 1]], [5, [-139189, -91471, -20563, -1457, 115, 23, 1]]], [[1], [0, 0, 1]]], [[[2, [-1, 1]], [6, [-31096, 18224, -3670, 121, 58, -13, 1]]], [[1], [0, 0, 1]]]]]

sage: print(generate_lines(db.artin_reps.lucky({'Baselabel':'2.66328117311920648.4t3.a'})))
2.66328117311920648.4t3.a.a:2:66328117311920648:[-95690836119520, -107985549736, -34145628, -1, 1]:[0, 0]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [-1]], [[1], [0], [-1]], [[1], [0], [1]], [[1], [-1]], [[1]]]:[[2, 17, 89, 617, 8681], [6, 7, 1, 7, 7]]:[[[1, 1, 1, 1], [1, 1, 2], [4]], [1, 4, 5]]:[[], []]:[[[2, 2]], [[[[2, [107996931612, 1]], [3, [11663337231837031108232, 215993863224, 1]]], [[1], [0, 0, 1]]]]]

sage: print(generate_lines(db.artin_reps.lucky({'Baselabel':'2.2e8_3e2.8t5.1'})))
2.2304.8t5.a.a:2:2304:[9, 0, -36, 0, 36, 0, -12, 0, 1]:[0, 0]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1]]]:[[2, 3], [6, 6]]:[[[1, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2]], [1, 2]]:[[], []]:[[[4, 4]], [[[[3, [-24470208, 0, 1]], [4, [-37324800, 0, 1]], [5, [-7312896, 0, 1]]], [[1], [8, 7, 6, 5, 4, 3, 2, 1]]]]]

sage: print(generate_lines(db.artin_reps.lucky({'Baselabel':'2.2e8_3e2.8t5.2'})))
2.2304.8t5.b.a:2:2304:[9, 0, 36, 0, 36, 0, 12, 0, 1]:[1, 1]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1]]]:[[2, 3], [6, 6]]:[[[1, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2]], [1, 2]]:[[], []]:[[[4, 4]], [[[[3, [-6015168, 0, 1]], [4, [-1672704, 0, 1]], [5, [-9331200, 0, 1]]], [[1], [8, 7, 6, 5, 4, 3, 2, 1]]]]]


# Running
# export_search({'Baselabel':{'$in':['2.47.5t2.a', '2.2304.8t5.a', '2.2304.8t5.b', '3.3003175.6t11.b','2.66328117311920648.4t3.a']}},'foo')
will generate the following lines
"""
2.47.5t2.a.a:2:47:[1, 0, -1, 2, -2, 1]:[0, 1]:5:[[[1], [-2], [1]], [[1], [0], [-1]], [[1], [0, 0, -1, -1], [1]], [[1], [1, 0, 1, 1], [1]], [[1], [-1]]]:[[2, 47], [3, 5]]:[[[1, 1, 1, 1, 1], [1, 2, 2]], [1, 2]]:[[[5]], [[[[4, [-47, 1]], [3, [47, 1]]], [2, 3, 4, 5, 1]]]]:[[], []]
2.47.5t2.a.b:2:47:[1, 0, -1, 2, -2, 1]:[0, 1]:5:[[[1], [-2], [1]], [[1], [0], [-1]], [[1], [1, 0, 1, 1], [1]], [[1], [0, 0, -1, -1], [1]], [[1], [-1]]]:[[2, 47], [3, 5]]:[[[1, 1, 1, 1, 1], [1, 2, 2]], [1, 2]]:[[[5]], [[[[4, [-47, 1]], [3, [47, 1]]], [2, 3, 4, 5, 1]]]]:[[], []]
2.2304.8t5.a.a:2:2304:[9, 0, -36, 0, 36, 0, -12, 0, 1]:[0, 0]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1]]]:[[2, 3], [6, 6]]:[[[1, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2]], [1, 2]]:[[], []]:[[[4, 4]], [[[[3, [-24470208, 0, 1]], [4, [-37324800, 0, 1]], [5, [-7312896, 0, 1]]], [[1], [8, 7, 6, 5, 4, 3, 2, 1]]]]]
2.2304.8t5.b.a:2:2304:[9, 0, 36, 0, 36, 0, 12, 0, 1]:[1, 1]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1], [0], [1]], [[1]]]:[[2, 3], [6, 6]]:[[[1, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2]], [1, 2]]:[[], []]:[[[4, 4]], [[[[3, [-6015168, 0, 1]], [4, [-1672704, 0, 1]], [5, [-9331200, 0, 1]]], [[1], [8, 7, 6, 5, 4, 3, 2, 1]]]]]
2.66328117311920648.4t3.a.a:2:66328117311920648:[-95690836119520, -107985549736, -34145628, -1, 1]:[0, 0]:1:[[[1], [-2], [1]], [[1], [2], [1]], [[1], [0], [-1]], [[1], [0], [-1]], [[1], [0], [1]], [[1], [-1]], [[1]]]:[[2, 17, 89, 617, 8681], [6, 7, 1, 7, 7]]:[[[1, 1, 1, 1], [1, 1, 2], [4]], [1, 4, 5]]:[[], []]:[[[2, 2]], [[[[2, [107996931612, 1]], [3, [11663337231837031108232, 215993863224, 1]]], [[1], [0, 0, 1]]]]]
3.3003175.6t11.b.a:3:3003175:[1, -1, 3, 3, 3, -1, 1]:[0, 0, 1]:1:[[[1], [-3], [3], [-1]], [[1], [3], [3], [1]], [[1], [-1], [-1], [1]], [[1], [1], [-1], [-1]], [[1], [1], [-1], [-1]], [[1], [-1], [-1], [1]], [[1], [0], [0], [-1]], [[1], [1], [1], [1]], [[1], [-1], [1], [-1]], [[1], [0], [0], [1]], [[1], [-1]], [[1], [0], [-1]]]:[[3, 5, 7, 131], [8, 11, 12, 11]]:[[[1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 2], [1, 1, 4], [2, 4], [3, 3], [6]], [1, 3, 8, 9, 7, 10]]:[[], []]:[[[1, 1, 2, 2], [2, 2, 2]], [[[[4, [-353, 22, 15, 1]], [5, [-139189, -91471, -20563, -1457, 115, 23, 1]]], [[1], [0, 0, 1]]], [[[2, [-1, 1]], [6, [-31096, 18224, -3670, 121, 58, -13, 1]]], [[1], [0, 0, 1]]]]]
"""
*/
