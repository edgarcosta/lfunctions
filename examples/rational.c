// Copyright Edgar Costa 2024
// See LICENSE file for license details.
//
//
/*
 * Computes the L-function data for a smalljac curve
 *
 *
 * Usage:
 *
 * Input file:
 * label:degree:conductor:weight:mus:euler_factors
 * For example an Elliptic curve curve:
 * 11a2:2:11:1:[0,1]:[[1,2,2],[1,1,3],[1,-1,5],[1,2,7],[1,-1,0],[1,-4,13],[1,2,17],[1,0,19],[1,1,23],[1,0,29],[1,-7,31],[1,-3,37],[1,8,41],[1,6,43],[1,-8,47],[1,6,53],[1,-5,59],[1,-12,61],[1,7,67],[1,3,71],[1,-4,73],[1,10,79]]
 * or Genus 2 curve, where we take advantage that we will only use a_n with n <= 138.84*sqrt(conductor)
 * 169.a.169.1:4:169:1:[0,0,1,1]:[[1,3,5,6,4],[1,2,1,6,9],[1,0,-7,0,25],[1,0,7,0,0],[1,0,11,0,0],[1,5,13,0,0],[1,-3,-8,0,0],[1,6,31,0,0],[1,-6,13,0,0],[1,3,-20,0,0],[1,0,-50,0,0],[1,-15,112,0,0],[1,9,68,0,0],[1,8,0,0,0],[1,0,0,0,0],[1,6,0,0,0],[1,-12,0,0,0],[1,1,0,0,0],[1,-6,0,0,0],[1,-6,0,0,0],[1,0,0,0,0],[1,-8,0,0,0],[1,0,0,0,0],[1,12,0,0,0],[1,-12,0,0,0],[1,-3,0,0,0],[1,20,0,0,0],[1,6,0,0,0],[1,0,0,0,0],[1,-15,0,0,0],[1,-2,0,0,0],[1,-36,0,0,0],[1,27,0,0,0],[1,-4,0,0,0],[1,33,0,0,0],[1,0,0,0,0],[1,26,0,0,0],[1,-36,0,0,0],[1,24,0,0,0],[1,6,0,0,0],[1,0,0,0,0],[1,-22,0,0,0],[1,18,0,0,0],[1,9,0,0,0],[1,-24,0,0,0],[1,-2,0,0,0],[1,10,0,0,0],[1,18,0,0,0],[1,-42,0,0,0],[1,0,0,0,0],[1,-12,0,0,0],[1,0,0,0,0],[1,3,0,0,0],[1,-18,0,0,0],[1,3,0,0,0],[1,12,0,0,0],[1,-6,0,0,0],[1,-36,0,0,0],[1,-7,0,0,0],[1,0,0,0,0],[1,4,0,0,0],[1,-9,0,0,0],[1,0,0,0,0],[1,60,0,0,0],[1,-20,0,0,0],[1,0,0,0,0],[1,48,0,0,0],[1,46,0,0,0],[1,-30,0,0,0],[1,-24,0,0,0],[1,-57,0,0,0],[1,0,0,0,0],[1,-22,0,0,0],[1,19,0,0,0],[1,42,0,0,0],[1,36,0,0,0],[1,18,0,0,0],[1,24,0,0,0],[1,3,0,0,0],[1,-27,0,0,0],[1,18,0,0,0],[1,0,0,0,0],[1,12,0,0,0],[1,-17,0,0,0],[1,-28,0,0,0],[1,24,0,0,0],[1,12,0,0,0],[1,-3,0,0,0],[1,-39,0,0,0],[1,0,0,0,0],[1,-24,0,0,0],[1,-42,0,0,0],[1,12,0,0,0],[1,12,0,0,0],[1,0,0,0,0],[1,36,0,0,0],[1,-33,0,0,0],[1,-18,0,0,0],[1,-16,0,0,0],[1,0,0,0,0],[1,44,0,0,0],[1,27,0,0,0],[1,0,0,0,0],[1,-42,0,0,0],[1,-80,0,0,0],[1,0,0,0,0],[1,36,0,0,0],[1,0,0,0,0],[1,-60,0,0,0],[1,25,0,0,0],[1,-34,0,0,0],[1,21,0,0,0],[1,39,0,0,0],[1,0,0,0,0],[1,84,0,0,0],[1,33,0,0,0],[1,-24,0,0,0],[1,18,0,0,0],[1,-30,0,0,0],[1,-12,0,0,0],[1,81,0,0,0],[1,-19,0,0,0],[1,12,0,0,0],[1,-42,0,0,0],[1,-24,0,0,0],[1,-36,0,0,0],[1,9,0,0,0],[1,-48,0,0,0],[1,64,0,0,0],[1,0,0,0,0],[1,36,0,0,0],[1,60,0,0,0],[1,-16,0,0,0],[1,-26,0,0,0],[1,-60,0,0,0],[1,-12,0,0,0],[1,60,0,0,0],[1,-66,0,0,0],[1,42,0,0,0],[1,33,0,0,0],[1,0,0,0,0],[1,72,0,0,0],[1,-4,0,0,0],[1,0,0,0,0],[1,25,0,0,0],[1,-78,0,0,0],[1,0,0,0,0],[1,-6,0,0,0],[1,28,0,0,0],[1,0,0,0,0],[1,-21,0,0,0],[1,27,0,0,0],[1,-20,0,0,0],[1,36,0,0,0],[1,28,0,0,0],[1,0,0,0,0],[1,22,0,0,0],[1,81,0,0,0],[1,-14,0,0,0],[1,0,0,0,0],[1,30,0,0,0],[1,-6,0,0,0],[1,0,0,0,0],[1,6,0,0,0],[1,-75,0,0,0],[1,0,0,0,0],[1,2,0,0,0],[1,17,0,0,0],[1,0,0,0,0],[1,-90,0,0,0],[1,0,0,0,0],[1,-63,0,0,0],[1,-30,0,0,0],[1,57,0,0,0],[1,-28,0,0,0],[1,39,0,0,0],[1,-30,0,0,0],[1,0,0,0,0],[1,-4,0,0,0],[1,-14,0,0,0],[1,0,0,0,0],[1,96,0,0,0],[1,-98,0,0,0],[1,0,0,0,0],[1,-60,0,0,0],[1,3,0,0,0],[1,44,0,0,0],[1,0,0,0,0],[1,39,0,0,0],[1,-36,0,0,0],[1,-31,0,0,0],[1,-30,0,0,0],[1,-64,0,0,0],[1,-48,0,0,0],[1,36,0,0,0],[1,-18,0,0,0],[1,0,0,0,0],[1,-29,0,0,0],[1,0,0,0,0],[1,60,0,0,0],[1,3,0,0,0],[1,2,0,0,0],[1,117,0,0,0],[1,-106,0,0,0],[1,78,0,0,0],[1,-63,0,0,0],[1,0,0,0,0],[1,24,0,0,0],[1,87,0,0,0],[1,4,0,0,0],[1,43,0,0,0],[1,132,0,0,0],[1,-62,0,0,0],[1,42,0,0,0],[1,-66,0,0,0],[1,0,0,0,0],[1,4,0,0,0],[1,-57,0,0,0],[1,-72,0,0,0],[1,0,0,0,0],[1,17,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[1,24,0,0,0],[1,48,0,0,0],[1,-10,0,0,0],[1,-57,0,0,0],[1,-24,0,0,0],[1,70,0,0,0],[1,0,0,0,0],[1,11,0,0,0],[1,20,0,0,0],[1,-18,0,0,0],[1,-114,0,0,0],[1,-4,0,0,0],[1,0,0,0,0],[1,0,0,0,0],[1,-57,0,0,0],[1,54,0,0,0],[1,60,0,0,0],[1,66,0,0,0],[1,-62,0,0,0],[1,-56,0,0,0],[1,-63,0,0,0],[1,-63,0,0,0],[1,72,0,0,0],[1,-18,0,0,0],[1,12,0,0,0],[1,-12,0,0,0],[1,54,0,0,0],[1,-60,0,0,0],[1,-36,0,0,0],[1,0,0,0,0],[1,-59,0,0,0],[1,-36,0,0,0],[1,84,0,0,0],[1,34,0,0,0],[1,-108,0,0,0],[1,-84,0,0,0],[1,27,0,0,0],[1,16,0,0,0],[1,-48,0,0,0],[1,0,0,0,0],[1,29,0,0,0],[1,123,0,0,0],[1,10,0,0,0],[1,72,0,0,0],[1,0,0,0,0],[1,-54,0,0,0],[1,-75,0,0,0],[1,-26,0,0,0],[1,0,0,0,0],[1,135,0,0,0],[1,50,0,0,0],[1,-14,0,0,0],[1,30,0,0,0],[1,-60,0,0,0],[1,0,0,0,0],[1,-15,0,0,0]]
 *
 *
 */

#include <acb_poly.h>
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>
#include <primesieve.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "glfunc.h"

typedef struct {
	// some string
	char* label;

	// L-function
  Lerror_t ecode;
  Lfunc_t L;
  double* mus;
  int weight;
  int degree;

	// the conductor of the L-function
  int64_t conductor;


	// euler factors
	// as N*(degree + 1) matrix
	int64_t* euler_factors;
  size_t size_euler_factors;
} Lfunc_rational;

typedef Lfunc_rational Lfunc_rational_t[1];

// provided by Drew Sutherland, 2024
static inline int atoii (long v[], int n, char *s)
{
  while (isspace(*s)) s++;
  char c = *s;
  if ( c == '[' || c == '(' || c == '{' ) { s++; while ( isspace(*s) ) s++; }
  if ( *s != '-' && *s != '+' && ! isdigit(*s) ) return 0;
  if ( !n ) return -1;
  int i = 0;
  while ( i < n ) {
    v[i] = atol(s); i++;
    if (*s=='-') s++;
    while ( isdigit(*s) ) s++;
    if ( *s=='.') { s++; while ( isdigit(*s) ) s++; }
    while ( isspace(*s) ) s++;
    if (*s != ',') break;
    s++; while (isspace(*s)) s++;
    if ( *s != '-' && *s != '+' && ! isdigit(*s) ) return -1;
  }
  if ( c == '[' && *s != ']' ) return -1;
  if ( c == '(' && *s != ')' ) return -1;
  if ( c == '{' && *s != '}' ) return -1;
  return i;
}

static inline int atoiii (long v[], int d[], int m, int n, char *s) // d has length m, v has length m*n, ith list has length d[i] stored at v+i*n
{
  char *t;
  int i;

  while (isspace(*s)) s++;
  if ( *s++ != '[' ) return 0;
  while (isspace(*s)) s++;
  if ( *s != '[' ) return 0;
  if (!m) return -1;
  d[0] = atoii(v,n,s);
  if ( d[0] < 0 ) return -1;
  for ( i = 1 ; i < m && (t = strchr(s,']')) ; i++ ) {
    s = t+1;
    while (isspace(*s)) s++;
    if ( *s == ']' ) break;
    if ( *s++ != ',' ) return -1;
    while (isspace(*s)) s++;
    if ( *s != '[' ) return -1;
    d[i] = atoii(v+i*n, n, s);
    if ( d[i] < 0 ) return -1;
  }
  return i;
}

static inline double atoff (double v[], int n, char *s)
{
  while (isspace(*s)) s++;
  char c = *s;
  if ( c == '[' || c == '(' || c == '{' ) { s++; while ( isspace(*s) ) s++; }
  if ( *s != '-' && *s != '+' && ! isdigit(*s) ) return 0;
  if ( !n ) return -1;
  int i = 0;
  while ( i < n ) {
    v[i] = atof(s); i++;
    if (*s=='-') s++;
    while ( isdigit(*s) ) s++;
    if ( *s=='.') { s++; while ( isdigit(*s) ) s++; }
    while ( isspace(*s) ) s++;
    if (*s != ',') break;
    s++; while (isspace(*s)) s++;
    if ( *s != '-' && *s != '+' && ! isdigit(*s) ) return -1;
  }
  if ( c == '[' && *s != ']' ) return -1;
  if ( c == '(' && *s != ')' ) return -1;
  if ( c == '{' && *s != '}' ) return -1;
  return i;
}

static inline size_t replace_char(char* str, const char d, const char n) {
  char *p;
  size_t count=0;
  p = str;
  // Count occurance of d in string
  while( (p=strchr(p, d)) != NULL ) {
    if(d != n) // compiler should be smart enough to figure out these are constants
      *p = n; // replace delimiter.
    p++; // Skip past our old delimiter
    count++;
  }
  return count;
}

// remove the nulls
static inline void replace_null(char *str, const char d, const size_t count) {
  char* p = str;
  for(size_t i = 0; i < count; ++i) {
    p = strchr(p, 0);
    *p = d; // replace null by delimiter.
    p++; // Skip past our old delimiter
  }
}

// this changes the string
static inline int split(char * str, char delim, char ***array, size_t *length ) {
  char *p;
  char **res;
  size_t count = replace_char(str, delim, 0); // replace delimiter with nulls

  // this splits the string in one more
  count++;

  // allocate dynamic array
  res = calloc( 1, count * sizeof(char *));
  if( !res ) return -1;

  p = str;
  for(size_t k=0; k <count; ++k ){
    if( *p ) res[k] = p;  // Copy start of string
    p = strchr(p, 0);    // Look for next null
    p++; // Start of next string
  }
  *array = res;
  *length = count;
  return 0;
}





void populate_local_factors(Lfunc_rational_t L) {
  size_t bound = Lfunc_nmax(L->L);
  // printf("bound = %ld\n", bound);
  size_t size;
  int64_t * primes = (int64_t*) primesieve_generate_primes(0, bound, &size, INT64_PRIMES);
  // FIXME check that we have enough euler factors
  // use Lfunc_reduce_nmax if euler_factors is too short
  assert(L->size_euler_factors >= size);

  size_t d = L->degree;
  acb_poly_t local_factor;
  acb_poly_init(local_factor);
  acb_poly_fit_length(local_factor, d + 1);


  for (size_t i = 0; i < size; ++i) {
    acb_poly_zero(local_factor);
    _acb_poly_set_length(local_factor, d + 1);
    for(size_t j = 0; j <= d; ++j) {
      acb_set_si(local_factor->coeffs + j, L->euler_factors[i*(d+1) + j]);
    }
    Lfunc_use_lpoly(L->L, primes[i], local_factor);
  }

  acb_poly_clear(local_factor);
}



void Lfunc_rational_clear(Lfunc_rational_t L) {
  if(L->label != NULL)
    free(L->label);
  if(L->L != NULL)
    Lfunc_clear(L->L);
  if(L->mus != NULL)
    free(L->mus);
  if(L->euler_factors != NULL)
    free(L->euler_factors);
}

void Lfunc_rational_init(Lfunc_rational_t L) {
  L->label = NULL;
  L->L = NULL;
  L->mus = NULL;
  L->euler_factors = NULL;
}

int Lfunc_rational_set_s(Lfunc_rational_t L, char *s) {

  char **tokens;
  size_t tokens_length;
  int status = 0;

  status = split(s, ':', &tokens, &tokens_length);
  // printf("tokens_length: %d\n", tokens_length);

  if(tokens_length != 6)
    status = -1;

  if(status != -1) {
    // read label
    size_t len = strlen(tokens[0]) + 1;
    L->label = (char *) malloc(len * sizeof(char));
    strncpy(L->label, tokens[0], len);
    // printf("label = %s\n", L->label);

    // printf("label = %s %d\n", tokens[1], strlen(tokens[1]));
    L->degree = atol(tokens[1]);
    // printf("degree = %d\n", L->degree);
    L->conductor = atol(tokens[2]);
    // printf("conductor = %ld\n", L->conductor);
    L->weight = atol(tokens[3]);
    // printf("weight = %d\n", L->weight);
    L->mus = (double *)malloc(L->degree*sizeof(double));
    status = atoff(L->mus, L->degree, tokens[4]);
  }


  if(status != -1) {
    // assuming a matrix as input
    // this just counts commas
    size_t entries = replace_char(tokens[5], ',', ',') + 1;
    assert(entries % (L->degree + 1) == 0);
    L->size_euler_factors = entries/(L->degree + 1);
    int* d;
    d = (int *)malloc(L->size_euler_factors * sizeof(int));
    for(size_t i = 0; i < L->size_euler_factors; ++i)
      d[i] = L->degree + 1;
    L->euler_factors = (int64_t *) malloc( L->size_euler_factors * (L->degree + 1) * sizeof(int64_t));
    status = atoiii(L->euler_factors, d, L->size_euler_factors, L->degree + 1, tokens[5]);
  }
  if(status != -1) {
    L->L = Lfunc_init(L->degree, L->conductor, L->weight*0.5, L->mus, &L->ecode);
    if(fatal_error(L->ecode)) {
      fprint_errors(stderr, L->ecode);
      status = -1;
    }
  }

  replace_null(s, ':', tokens_length - 1);
  return status;
}




int main(int argc, char** argv) {
  assert(argc == 3);
  printf("Input: %s\n", argv[1]);
  printf("Output: %s\n", argv[2]);

  FILE* input = fopen(argv[1], "r");
  if (input == NULL) {
    printf("Could not open file %s.\n", argv[1]);
    return 1;
  }
  FILE* output = fopen(argv[2], "w");
  if (output == NULL) {
    printf("Could not open file %s.\n", argv[2]);
    return 1;
  }


  char *line = NULL;
  size_t len = 0;
  while (getline(&line, &len, input) != -1) {
    Lfunc_rational_t L;
    Lfunc_rational_init(L);
    Lfunc_rational_set_s(L, line);
    printf("label = %s\n", L->label);
    printf("degree = %d conductor = %" PRId64 " weight = %d mus = [ ", L->degree, L->conductor, L->weight);
    for(int i=0; i < L->degree; ++i)
      printf("%.2f ", L->mus[i]);
    printf("]\n");


    populate_local_factors(L);
    if(fatal_error(L->ecode)) {
      fprint_errors(stderr, L->ecode);
      return -1;
    }
    L->ecode|=Lfunc_compute(L->L);
    if(fatal_error(L->ecode)) {
      fprint_errors(stderr, L->ecode);
      return -1;
    }
    printf("Error code = %" PRIu64 "\n", L->ecode);
    printf("Rank = %" PRIu64 "\n", Lfunc_rank(L->L));
    printf("Epsilon = ");acb_printd(Lfunc_epsilon(L->L),20);printf("\n");
    printf("Leading Taylor coeff = ");arb_printd(Lfunc_Taylor(L->L), 20);printf("\n");
    printf("First zero = ");arb_printd(Lfunc_zeros(L->L, 0), 20);printf("\n");
    printf("\n\n");

    // TODO write output to output
    Lfunc_rational_clear(L);
  }
  fclose(input);
  fclose(output);
}
