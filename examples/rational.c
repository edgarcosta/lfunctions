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

#include <flint/acb_poly.h>
#include <flint/fmpz_poly.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
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

static bool checked_mul_size(size_t a, size_t b, size_t *out) {
  if(a != 0 && b > SIZE_MAX/a)
    return false;
  *out = a*b;
  return true;
}

static bool parse_int64_field(const char *s, int64_t *out) {
  char *end = NULL;
  errno = 0;
  while(isspace((unsigned char)*s))
    s++;
  intmax_t v = strtoimax(s, &end, 10);
  if(errno != 0 || end == s || v < INT64_MIN || v > INT64_MAX)
    return false;
  while(isspace((unsigned char)*end))
    end++;
  if(*end != '\0')
    return false;
  *out = (int64_t)v;
  return true;
}

static bool parse_int64_token(char **s, int64_t *out) {
  char *p = *s;
  char *end = NULL;
  errno = 0;
  intmax_t v = strtoimax(p, &end, 10);
  if(errno != 0 || end == p || v < INT64_MIN || v > INT64_MAX)
    return false;
  *out = (int64_t)v;
  *s = end;
  return true;
}

static bool parse_double_token(char **s, double *out) {
  char *p = *s;
  char *end = NULL;
  errno = 0;
  double v = strtod(p, &end);
  if(errno != 0 || end == p || !isfinite(v))
    return false;
  *out = v;
  *s = end;
  return true;
}

// provided by Drew Sutherland, 2024
static inline int atoii (int64_t v[], int n, char *s)
{
  while (isspace((unsigned char)*s)) s++;
  char c = *s;
  if ( c == '[' || c == '(' || c == '{' ) { s++; while ( isspace((unsigned char)*s) ) s++; }
  if ( *s != '-' && *s != '+' && ! isdigit((unsigned char)*s) ) return 0;
  if ( !n ) return -1;
  int i = 0;
  while ( i < n ) {
    if(!parse_int64_token(&s, v+i))
      return -1;
    i++;
    while ( isspace((unsigned char)*s) ) s++;
    if (*s != ',') break;
    s++; while (isspace((unsigned char)*s)) s++;
    if ( *s != '-' && *s != '+' && ! isdigit((unsigned char)*s) ) return -1;
  }
  if ( c == '[' && *s != ']' ) return -1;
  if ( c == '(' && *s != ')' ) return -1;
  if ( c == '{' && *s != '}' ) return -1;
  return i;
}

static inline int atoiii (int64_t v[], int d[], int m, int n, char *s) // d has length m, v has length m*n, ith list has length d[i] stored at v+i*n
{
  char *t;
  int i;

  while (isspace((unsigned char)*s)) s++;
  if ( *s++ != '[' ) return 0;
  while (isspace((unsigned char)*s)) s++;
  if ( *s != '[' ) return 0;
  if (!m) return -1;
  d[0] = atoii(v,n,s);
  if ( d[0] < 0 ) return -1;
  for ( i = 1 ; i < m && (t = strchr(s,']')) ; i++ ) {
    s = t+1;
    while (isspace((unsigned char)*s)) s++;
    if ( *s == ']' ) break;
    if ( *s++ != ',' ) return -1;
    while (isspace((unsigned char)*s)) s++;
    if ( *s != '[' ) return -1;
    d[i] = atoii(v+i*n, n, s);
    if ( d[i] < 0 ) return -1;
  }
  return i;
}

static inline int atoff (double v[], int n, char *s)
{
  while (isspace((unsigned char)*s)) s++;
  char c = *s;
  if ( c == '[' || c == '(' || c == '{' ) { s++; while ( isspace((unsigned char)*s) ) s++; }
  if ( *s != '-' && *s != '+' && ! isdigit((unsigned char)*s) ) return 0;
  if ( !n ) return -1;
  int i = 0;
  while ( i < n ) {
    if(!parse_double_token(&s, v+i))
      return -1;
    i++;
    while ( isspace((unsigned char)*s) ) s++;
    if (*s != ',') break;
    s++; while (isspace((unsigned char)*s)) s++;
    if ( *s != '-' && *s != '+' && ! isdigit((unsigned char)*s) ) return -1;
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
  if( !res ) {
    replace_null(str, delim, count - 1);
    return -1;
  }

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
  size_t n = L->size_euler_factors;
  size_t d = L->degree;
  if (n == 0) {
    fprintf(stderr, "No Euler factors parsed for %s.\n", L->label ? L->label : "(unknown)");
    L->ecode |= ERR_NO_DATA;
    return;
  }

  // Pack the parsed integer Euler factors into one contiguous fmpz_poly array,
  // factor i (= the i-th prime, 2, 3, 5, ...) being 1 + c_1 T + ... + c_d T^d.
  fmpz_poly_struct *factors = (fmpz_poly_struct *) malloc(n * sizeof(fmpz_poly_struct));
  if(!factors) {
    L->ecode |= ERR_OOM;
    return;
  }
  for (size_t i = 0; i < n; ++i) {
    fmpz_poly_init(factors + i);
    for (size_t j = 0; j <= d; ++j)
      fmpz_poly_set_coeff_si(factors + i, j, L->euler_factors[i*(d+1) + j]);
  }

  // One array call replaces the manual prime sieve and the per-prime push: the
  // library sieves the primes itself and consumes the k-th factor for the k-th
  // prime. If the array is shorter than the library needs it reduces nmax and
  // warns (ERR_INSUFF_EULER) rather than the old hard assert, so a truncated
  // factor list (as in the genus-2 example above) just works.
  L->ecode |= Lfunc_use_lpolys_fmpz(L->L, factors, n);

  for (size_t i = 0; i < n; ++i)
    fmpz_poly_clear(factors + i);
  free(factors);
}



void Lfunc_rational_init(Lfunc_rational_t L) {
  L->label = NULL;
  L->ecode = ERR_SUCCESS;
  L->L = NULL;
  L->mus = NULL;
  L->weight = 0;
  L->degree = 0;
  L->conductor = 0;
  L->euler_factors = NULL;
  L->size_euler_factors = 0;
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
  Lfunc_rational_init(L);
}

int Lfunc_rational_set_s(Lfunc_rational_t L, char *s) {

  char **tokens = NULL;
  size_t tokens_length = 0;
  size_t label_len = 0;
  int status = 0;

  status = split(s, ':', &tokens, &tokens_length);
  // printf("tokens_length: %d\n", tokens_length);

  if(status != 0)
    return -1;

  if(tokens_length != 6)
    status = -1;

  if(status != -1) {
    for(size_t i = 0; i < tokens_length; ++i) {
      if(tokens[i] == NULL) {
        status = -1;
        break;
      }
    }
  }

  if(status != -1) {
    // read label
    label_len = strlen(tokens[0]) + 1;
    L->label = (char *) malloc(label_len * sizeof(char));
    if(!L->label) {
      L->ecode |= ERR_OOM;
      status = -1;
    }
  }

  if(status != -1) {
    strncpy(L->label, tokens[0], label_len);
    // printf("label = %s\n", L->label);

    int64_t degree64 = 0, conductor64 = 0, weight64 = 0;
    if(!parse_int64_field(tokens[1], &degree64) ||
       !parse_int64_field(tokens[2], &conductor64) ||
       !parse_int64_field(tokens[3], &weight64) ||
       degree64 <= 0 || degree64 > MAX_DEGREE ||
       conductor64 <= 0 ||
       weight64 < INT_MIN || weight64 > INT_MAX) {
      status = -1;
    }
    else {
      L->degree = (int)degree64;
      L->conductor = conductor64;
      L->weight = (int)weight64;
    }
  }

  if(status != -1) {
    size_t bytes = 0;
    if(!checked_mul_size((size_t)L->degree, sizeof(double), &bytes)) {
      status = -1;
    }
    else {
      L->mus = (double *)malloc(bytes);
    }
    if(!L->mus) {
      L->ecode |= ERR_OOM;
      status = -1;
    }
  }

  if(status != -1) {
    status = atoff(L->mus, L->degree, tokens[4]);
    if(status != L->degree)
      status = -1;
  }


  if(status != -1) {
    // assuming a matrix as input
    // this just counts commas
    size_t entries = replace_char(tokens[5], ',', ',') + 1;
    size_t width = (size_t)L->degree + 1;
    if(entries % width != 0) {
      status = -1;
    }
    else {
      L->size_euler_factors = entries/width;
      if(L->size_euler_factors == 0 || L->size_euler_factors > (size_t)INT_MAX)
        status = -1;
    }
  }

  if(status != -1) {
    int* d = NULL;
    size_t width = (size_t)L->degree + 1;
    size_t entries = 0, bytes = 0;
    if(!checked_mul_size(L->size_euler_factors, sizeof(int), &bytes)) {
      status = -1;
    }
    if(status != -1)
      d = (int *)malloc(bytes);
    if(!d) {
      L->ecode |= ERR_OOM;
      status = -1;
    }
    if(status != -1) {
      for(size_t i = 0; i < L->size_euler_factors; ++i)
        d[i] = (int)width;
      if(!checked_mul_size(L->size_euler_factors, width, &entries) ||
         !checked_mul_size(entries, sizeof(int64_t), &bytes)) {
        status = -1;
      }
    }
    if(status != -1) {
      L->euler_factors = (int64_t *) malloc(bytes);
      if(!L->euler_factors) {
        L->ecode |= ERR_OOM;
        status = -1;
      }
    }
    if(status != -1) {
      status = atoiii(L->euler_factors, d, (int)L->size_euler_factors, (int)width, tokens[5]);
      if(status <= 0 || (size_t)status != L->size_euler_factors) {
        L->size_euler_factors = status > 0 ? (size_t)status : 0;
        status = -1;
      }
    }
    if(status != -1) {
      for(size_t i = 0; i < L->size_euler_factors; ++i) {
        if(d[i] != L->degree + 1) {
          status = -1;
          break;
        }
      }
    }
    free(d);
  }
  if(status != -1) {
    L->L = Lfunc_init(L->degree, L->conductor, L->weight*0.5, L->mus, &L->ecode);
    if(fatal_error(L->ecode)) {
      status = -1;
    }
  }

  replace_null(s, ':', tokens_length - 1);
  free(tokens);
  if(status == -1)
    L->size_euler_factors = 0;
  return status == -1 ? -1 : 0;
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
    fclose(input);
    return 1;
  }


  char *line = NULL;
  size_t len = 0;
  size_t line_no = 0;
  int rc = 0;
  while (getline(&line, &len, input) != -1) {
    line_no++;
    Lfunc_rational_t L;
    Lfunc_rational_init(L);
    if(Lfunc_rational_set_s(L, line) != 0) {
      fprintf(stderr, "Could not parse input line %zu.\n", line_no);
      if(fatal_error(L->ecode))
        fprint_errors(stderr, L->ecode);
      Lfunc_rational_clear(L);
      rc = 1;
      break;
    }
    printf("label = %s\n", L->label);
    printf("degree = %d conductor = %" PRId64 " weight = %d mus = [ ", L->degree, L->conductor, L->weight);
    for(int i=0; i < L->degree; ++i)
      printf("%.2f ", L->mus[i]);
    printf("]\n");


    populate_local_factors(L);
    if(fatal_error(L->ecode)) {
      fprint_errors(stderr, L->ecode);
      Lfunc_rational_clear(L);
      rc = 1;
      break;
    }
    L->ecode|=Lfunc_compute(L->L);
    if(fatal_error(L->ecode)) {
      fprint_errors(stderr, L->ecode);
      Lfunc_rational_clear(L);
      rc = 1;
      break;
    }
    if(L->ecode != ERR_SUCCESS) {
      fprintf(stderr, "Warnings for %s:\n", L->label);
      fprint_errors(stderr, L->ecode);
    }
    printf("Rank = %" PRIu64 "\n", Lfunc_rank(L->L));
    printf("Epsilon = ");acb_printd(Lfunc_sign(L->L),20);printf("\n");
    printf("Leading Taylor coeff = ");arb_printd(Lfunc_Taylor(L->L), 20);printf("\n");
    printf("First zero = ");arb_printd(Lfunc_zeros(L->L, 0), 20);printf("\n");

    // TODO write output to output
    Lfunc_rational_clear(L);
  }
  free(line);
  fclose(input);
  fclose(output);
  return rc;
}
