// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Some of the common tools used in the examples
 */


#define __STDC_FORMAT_MACROS
#include <cassert>
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
#include <acb.h>


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
using std::multimap;
using std::ostream;
using std::ofstream;
using std::pair;
using std::runtime_error;
using std::set;
using std::size_t;
using std::string;
using namespace std::string_literals;
using std::stringstream;
using std::vector;

using flint::fmpzxx;


/******************************************************************************
 * Generic handy tools
 *  - my_exception: an exception class that throws line number
 *  - assert_print: a better assert
 *  - print: lazy print
 *****************************************************************************/

// exception that throws line number
class my_exception : public runtime_error {
    string msg;
public:
    my_exception(const std::string &arg, const char *file, int line) :
    std::runtime_error(arg) {
        std::ostringstream o;
        o << file << ":" << line << ": " << arg;
        msg = o.str();
    }
    ~my_exception() throw() {}
    const char *what() const throw() {
        return msg.c_str();
    }
};
#define throw_line(arg) throw my_exception(arg, __FILE__, __LINE__);

// a better assert
#ifndef NDEBUG
#define assert_print(left, operator, right) \
{ \
    if( !( (left) operator (right) ) ) \
    { \
        cerr << "ASSERT FAILED: " << #left << " " << #operator << " " << #right << " @ " << __FILE__ << ":" << __LINE__  << endl; \
        cerr << #left << " = " << (left) << "; " << #right << " = " << (right) << endl; \
        abort(); \
    } \
}
#else
#define assert_print(condition, statement) ((void)0)
#endif

// lazy print
#define print(var) { cout << #var << " = " << (var) << endl;}

/******************************************************************************
 * in/out operators for various (templated) classes
 *  - operator >> for fmpzxx
 *  - operators << and >> for vector<T>
 *  - operators << and >> for map<T, R>
 *****************************************************************************/

//operator >> for fmpzxx
istream & operator>>(istream& s, fmpzxx& output) {
  stringstream buffer;
  long c;
  if (!s)
    throw_line("bad fmpzxx input"s);

  c = s.peek();
  while (iswspace(c)) {
    s.get();
    c = s.peek();
  }
  if (c != '-' and not iswdigit(c))
    throw_line("bad fmpzxx input"s);

  //puts the first character
  buffer.put(s.get());
  c = s.peek();
  while (iswdigit(c)) {
    buffer.put(s.get());
    c = s.peek();
  }
  mpz_t tmp;
  mpz_init(tmp);
  gmp_sscanf(buffer.str().c_str(),"%Zd", tmp);
  fmpz_set_mpz(output._fmpz(), tmp);
  mpz_clear(tmp);
  return s;
}


// outputs [a, b, c, ..]
template<class T>
ostream& operator<<(ostream& s, const vector<T>& a) {
  size_t n = a.size();
  s <<"[";
  for(size_t i = 0; i < n; ++i) {
    s << a[i];
    if(i < n - 1) s<<", ";
  }
  s << "]";
  return s;
}

// reads [a, b, c, ..]
template<class T>
istream & operator>>(istream& s, vector<T>& output) {
  vector<T> ibuf(0);
  long c;
  if (!s)
    throw_line("bad vector input"s);

  c = s.peek();
  while (iswspace(c)) {
    s.get();
    c = s.peek();
  }
  if (c != '[')
    throw_line("bad vector input"s);

  s.get();
  c = s.peek();
  while (iswspace(c)) {
    s.get();
    c = s.peek();
  }
  while (c != ']' and c != EOF) {
    T tmp;
    if (!(s >> tmp))
      throw_line("bad vector input"s);
    ibuf.push_back(tmp);
    c = s.peek();
    while (iswspace(c) or c == ',') {
      s.get();
      c = s.peek();
    }
  }

  if (c == EOF)
    throw_line("bad vector input"s);
  s.get();
  output = ibuf;
  return s;
}



// operator << for pair<T, R>
// outputs [a, b]
template<class T, class R>
ostream& operator<<(ostream& s, const pair<T, R>& a) {
  s <<"["<<a.first<<", "<<a.second<<"]";
  return s;
}

// reads [a, b]
template<class T, class R>
istream & operator>>(istream& s, pair<T, R>& output) {
  long c;
  if (!s)
    throw_line("bad pair input"s);

  c = s.peek();
  while (iswspace(c)) {
    s.get();
    c = s.peek();
  }
  if (c != '[')
    throw_line("bad pair input"s);

  s.get();
  c = s.peek();
  while (iswspace(c)) {
    s.get();
    c = s.peek();
  }
  if(!(s>>output.first))
    throw_line("bad first element"s);

  c = s.peek();
  while (iswspace(c) or c == ',') {
    s.get();
    c = s.peek();
  }

  if(!(s>>output.second))
    throw_line("bad second element"s);

  while (c != ']' and c != EOF) {
    c = s.peek();
  }

  if (c == EOF)
    throw_line("bad pair input"s);
  s.get();

  return s;
}

/* operator >> for (multi)map<T, R>
 * reads a vector with keys and vector with values
 * and returns a map
 * [[k1, k2, ...], [v1, v2, ...]]
 * correponds to
 * k1 -> v1
 * k2 -> v2
 * ...
 */
template<typename T, typename R, typename Compare, typename Allocator>
istream & operator>>(istream& s, map<T, R, Compare, Allocator>& a)
{
  pair<vector<T>, vector<R>> kv;
  s >> kv;
  vector<T> &keys = kv.first;
  vector<R> &values = kv.second;

  assert_print(keys.size(), ==, values.size());
  map<T, R, Compare, Allocator> ibuf;
  for(size_t i = 0; i < values.size(); ++i)
    ibuf.emplace(make_pair(keys[i], values[i]));
  a = ibuf;
  return s;
}
template<typename T, typename R, typename Compare, typename Allocator>
istream & operator>>(istream& s, multimap<T, R, Compare, Allocator>& a)
{
  pair<vector<T>, vector<R>> kv;
  s >> kv;
  vector<T> &keys = kv.first;
  vector<R> &values = kv.second;

  assert_print(keys.size(), ==, values.size());
  multimap<T, R, Compare, Allocator> ibuf;
  for(size_t i = 0; i < values.size(); ++i)
    ibuf.emplace(make_pair(keys[i], values[i]));
  a = ibuf;
  return s;
}




/* outputs a python dictionary
 * {k1: v1
 *  k2: v2
 *  ...}
 */
template<class T, class R, class Compare>
ostream & operator<<(ostream& s, const map<T, R, Compare>& a)
{
    s << "{";
    int64_t i, n;
    i = 0;
    n = a.size();
    for(auto it = a.cbegin() ; a.cend() != it ; ++it) {
        s << it->first;
        s << ": ";
        s << it->second;
        if(i < n - 1)
            s <<",\n ";
        /*else
            break;*/
        ++i;
    }
    s << "}";
    return s;
}

// retuns the interval [a*2^e, b*2^e] as [a, b, e]
ostream& operator<<(ostream& s, const arb_t x) {
  vector<fmpzxx> tmp(3);
  fmpzxx &a = tmp[0];
  fmpzxx &b = tmp[1];
  fmpzxx &e = tmp[2];
  arb_get_interval_fmpz_2exp(a._fmpz(), b._fmpz(), e._fmpz(), x);
  s << tmp;
  return s;
}

// returns the rectangle [ar*2^er, br*2^er] + [ai*2^ei, bi*2^ei]*I
// as [[ar, br, er], [ai, bi, ei]]
ostream& operator<<(ostream &s, const acb_t z) {
  s << "[" << acb_realref(z) <<", "<< acb_imagref(z) << "]";
  return s;
}


// << operator for *Lplot_t
ostream& operator<<(ostream& s, const Lplot_t *Lpp) {
  s << Lpp->spacing << ":";
  vector<double> plot(Lpp->points,Lpp->points + Lpp->n_points);
  s << plot;
  return s;
}


/*
 * GCD and LCM for int64_t using n_gcd_full from flint
 */
inline int64_t gcd(const int64_t a, const int64_t b) {
  return int64_t(n_gcd_full(
        a >= 0 ? mp_limb_t(a) : mp_limb_t(-a),
        b >= 0 ? mp_limb_t(b) : mp_limb_t(-b)));
}


int64_t lcm(int64_t a, int64_t b) {
  return a*b/gcd(a, b);
}
int64_t lcm(const vector<size_t>& v) {
  int64_t res = 1;
  for(const auto &elt : v)
    res = lcm(res, int64_t(elt));
  return res;
}
