// Copyright Edgar Costa 2019
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
 * For example:
 *
 */

#include <stdio.h>

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
  replace = d == n;
  // Count occurance of d in string
  while( (p=strchr(p, d)) != NULL ) {
    if(replace)
      *p = n; // replace delimiter.
    p++; // Skip past our old delimiter
    count++;
  }
  return count;
}

// remove the nulls
static inline replace_null(char *str, const char d, const size_t count) {
  char* p = str;
  for(size_t i = 0; i < count; ++i) {
    p=strchr(p, 0);
    *p = n; // replace delimiter.
    p++; // Skip past our old delimiter
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





long populate_local_factors(Lfunc_rational_t L) {
  size_t bound = Lfunc_nmax(L->L);
  size_t size;
  int64_t * primes = (int64_t*) primesieve_generate_primes(0, bound, &size, INT64_PRIMES);
  // FIXME check that we have enough euler factors
  // use Lfunc_reduce_nmax if euler_factors is too short
  assert size_euler_factors
  size = max(size, size_euler_factors);

  size_t d = L->degree;
  acb_poly_t local_factor;
  acb_poly_init(local_factor);
  acb_poly_fit_length(local_factor, d + 1);



  for (size_t i = 0; i < size; ++i) {
    int64_t p = primes[i];
    int64_t pk = 1;
    acb_poly_zero(local_factor);
    for(size_t j=0; j <= d; ++j) {
      acb_set_d(local_factor->coeffs + i, L->euler_factors + i*(d+1) + j);
    }
    Lfunc_use_lpoly(L->L, p, local_factor);
  }

  primesieve_free_iterator(&it);
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

  if(tokens_length != 6)
    status = -1;

  if(status != -1) {
    // read label
    L->label = (char *) malloc((strlen(tokens[0]) + 1) * sizeof(char));
    strncpy(L->label, input, len);


    L->degree = atol(tokens[1]);
    L->conductor = atol(tokens[2]);
    L->weight = atol(tokens[3]);
    status = atoff(L->mus, L->degree, tokens[4]);


  }


  if(status != -1) {
    // assuming a matrix as input
    size_t entries = replace_char(tokens[5], ',', ',');
    size_t L->size_euler_factors = entries/(d + 1);
    int* d;
    d = (int *)malloc(L->size_euler_factors * sizeof(int));
    for(size_t i = 0; i < L->size_euler_factors; ++i)
      d[i] = L->degree + 1;
    status = aoiii(L->euler_factors, d, L->size_euler_factors, L->degree + 1, input)
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
  while ((ssize_t read = getline(&line, &len, file)) != -1) {
    printf("READ: %s", line);
    Lfunc_rational L;
    Lfunc_rational_init(L);
    Lfunc_rational_set_s(L, s);


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
    printf("Rank = %" PRIu64 "\n",Lfunc_rank(L));
    printf("Epsilon = ");acb_printd(Lfunc_epsilon(L),20);printf("\n");
    printf("Leading Taylor coeff = ");arb_printd(Lfunc_Taylor(L), 20);printf("\n");
    printf("First zero = ");arb_printd(Lfunc_zeros(L, 0), 20);printf("\n");

    Lfunc_rational_clear(L);
  }


  fclose(input);
  fclose(output);
}
