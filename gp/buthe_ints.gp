\p400
allocatemem(2^30);
MAX_MU=100;
MAX_MU2=200;
MAX_R=9;
h=4;
bint(r,mu)=2*intnum(t=0,100,exp(-(0.5+mu)*t)/(1-exp(-2*t))*(64/(r*Pi)-sin(64*t/r)/(Pi*t*cosh(2*t))));

trunc(x)=round(x*2^20); /* +/- 1 added in buthe.h */

printf("{");
for(r=2,2,for(kk=0,2*MAX_MU2,k=kk/2.0;b=64.0/r;S=bint(r,k);printf("%d,",trunc(S)));printf("\n");)

for(r=3,MAX_R,for(kk=0,2*MAX_MU,k=kk/2.0;b=64.0/r;S=bint(r,k);printf("%d,",trunc(S)));for(kk=2*MAX_MU+1,2*MAX_MU2,printf("0,"));printf("\n");)
printf("}\n");
quit;
