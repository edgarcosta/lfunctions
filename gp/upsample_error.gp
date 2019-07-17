/p 200
check_M(M)=M>0.5;

check_H(H,h,A,mus,r)=
{
  s=sum(d=1,r,mus[d]);
  s+=r/2.0;
  s/=Pi;
  s=sqrt(s*A*h/2.0);
  return(H>s);
}

check_l(M,H,A,l)=
{
  tmp=M*A;
  if(H<tmp,tmp=H);
  return(l<Pi*tmp);
}

E(s,t,h,mus,r,N,T)=
{
  res=1.0;
  for(j=1,r,res*=exp((s+mus[j]-0.5)/2.0*log(sqrt((s+mus[j]-0.5)*(s+mus[j]-0.5)+(T+t)*(T+t))/(2.0*Pi))));
  res*=exp(-Pi*t*t/(h*h));
  res*=exp(s/2*log(N));
  res*=exp(r/2.0*log(2.0));
  res*=exp(r*log(zeta(s+0.5)));
  return(res);
}
upsample_error(M,H,h,A,mus,r,N,T,imz,l)=
{
  if(!check_M(M),return(0.0));
  if(!check_H(H,h,A,mus,r),return(0.0));
  if(!check_l(M,H,A,l),return(0.0));
  r0=0;
  for(d=1,r,if(mus[d]<0.5,r0++));
  delta=0.0;
  for(j=1,r,delta+=(M+mus[j]-0.5)/((M+mus[j]-0.5)*(M+mus[j]-0.5)+T*T));
  delta*=h*h;
  delta/=4*Pi;
  if(delta>=1.0,return(0.0));
  printf("delta=%10.8e\n",delta);
  E0=E(M,0,h,mus,r,N,T);
  E1=E(1,H/A,h,mus,r,N,T);
  E2=E(1,(H+1)/A,h,mus,r,N,T);
print(E0);print(E1);print(E2);
  printf("E0=%10.8e E1=%10.8e E2=%10.8e\n",E0,E1,E2);
  return(2*cosh(Pi*A*imz)*exp(l*log(Pi*A))*(A*h*exp(Pi/(h*h)*(M*M+(delta*T)*(delta*T)/1-delta)-Pi*M*A)/((Pi*M*A-l)*sqrt(1-delta))*E0+exp(0.75*r0*log(3.0))*E1*E1/((Pi*H-l)*(E1-E2))));
}

stride=4;
r=2;
mus=[0,1];
T=512.0/r/8;
  
A=2^16/2^9*r/stride;
  
h=sqrt(1.0/A)*1.001;
H=ceil(A*A*h*h/2.0);
M=H/A;
im=-0.05;
N=87;
err=upsample_error(M,H,h,A,mus,r,N,T,im,0);
printf("h=%10.8e M=%10.8e H=%10.8e\n",h,M,H);
printf("upsample error = %10.8e\n",err);
while(-log(err)/log(2)<100,h*=1.01;H=ceil(A*A*h*h/2.0);M=H/A;err=upsample_error(M,H,h,A,mus,r,N,T,im,0);printf("h=%10.8e M=%10.8e H=%10.8e\n",h,M,H);printf("upsample error = %10.8e\n",err));
