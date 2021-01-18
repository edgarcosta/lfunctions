#ifndef _M61MOD_INCLUDE_
#define _M61MOD_INCLUDE_

#define M61LOW		0x3FFFFFFFUL
#define M61LOWBITS	30
#define M61HIGH		0x7FFFFFFFUL
#define M61HIGHBITS	31
#define M61UI		0x1FFFFFFFFFFFFFFFUL		// this should raise a compile-time error on 32-bit machines
#define M61SI		0x1FFFFFFFFFFFFFFFL

#ifdef __cplusplus
extern "C"{
#endif

static inline long m61mod_ui (unsigned long n)
	{ while ( n > M61UI ) n = (n&M61UI)+(n>>61); return n == M61UI ? 0 : n; }
static inline long m61mod_si (long n)
	{ if ( n >= 0 ) return m61mod_ui(n); n = m61mod_ui(-n);  return n ? (long)M61SI-n : 0; }


// m61add and m61mul take inputs and return output in [0,2^61-1]
static inline long m61add (long a, long b)
	{ long c = a + b;  if ( c >= M61SI ) c-= M61SI; return c; }
static inline long m61mul (long a, long b)
{
	register unsigned long a0 = ((unsigned long)a)&M61LOW;
	register unsigned long a1 = ((unsigned long)a)>>M61LOWBITS;
	register unsigned long b0 = ((unsigned long)b)&M61HIGH;
	register unsigned long b1 = ((unsigned long)b)>>M61HIGHBITS;
	register unsigned long c = a0*b0 + a1*b1;
	register unsigned long d = ((a0*b1)<<1) + a1*b0;
	c += ((d<<M61LOWBITS)&M61UI) + (d>>M61HIGHBITS);
	return m61mod_ui (c);
}
#ifdef __cplusplus
}
#endif
#endif
