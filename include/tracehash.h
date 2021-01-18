#ifndef _TRACEHASH_INCLUDE_
#define _TRACEHASH_INCLUDE_

#define APHASH_MINP             4099    // next_prime(2^12)
#define APHASH_MINPI            565     // pi(APHASH_MINP)
#define APHASH_MAXP             8191    // prev_prime(2^13)
#define APHASH_MAXPI            1028    // pi(APHASH_MAXP)
#define APHASH_PRIMES           (APHASH_MAXPI-APHASH_MINPI+1)

#ifdef __cplusplus
extern "C"{
#endif

// compute trace hash given list of ap's for p in [APHASH_MINP,APHASH_MAXP] (so ap[0] should be a_4099)
long aphash_list (long ap[APHASH_PRIMES]);
// the ap should be either long or reduced 2^61-1
long aphash (long (*ap)(long p, void *arg), void *arg);

#ifdef __cplusplus
}
#endif
#endif
