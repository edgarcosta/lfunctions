/* No longer compiled into the library: src/buthe_winf.c now computes this
   integral rigorously on the fly (see buthe_winf_integral), so the captured
   buthe_ints.out table is no longer #included. Kept as a dev cross-check for
   the C routine.

   This Pari code computes an approximation to the integral in W_\infty for
   Buthe's zero-checking algorithm.
*/
\p400 /* lots of bits of precision please */
allocatemem(2^30); /* lots of memory too */
MAX_MU=100;
MAX_MU2=200;
MAX_R=9;
h=8;

/* compute 2\int\limits_0^{100}\frac{\exp(-(1/2-\mu_j)t)}{1-\exp(-2t)}(f(0)-f(t)) dt
   assumes 1) Pari can numerically integrate this to better than 20 bits
           2) The tail from 100->infty is negligible
*/

bint(b,mu)=2*intnum(t=0,100,exp(-(0.5+mu)*t)/(1-exp(-2*t))*(b/Pi-sin(b*t)/(Pi*t*cosh(h*t/2))));

trunc(x)=round(x*2^20); /* +/- 1 added in buthe.h */

printf("{");
for(r=2,2,for(kk=0,2*MAX_MU2,k=kk/2.0;b=64.0/r;S=bint(b,k);printf("%d,",trunc(S)));printf("\n");)

for(r=3,MAX_R,for(kk=0,2*MAX_MU,k=kk/2.0;b=64.0/r;S=bint(b,k);printf("%d,",trunc(S)));for(kk=2*MAX_MU+1,2*MAX_MU2,printf("0,"));printf("\n");)
printf("}\n");

printf("#define BUTHE_H (%d)\n",h);
quit;
