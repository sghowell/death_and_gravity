## Conditioned cloud marks and explicit moments

Use only S329's FIXED BORN positive D4 leading cloud of index
0<=a<=1, and total cut R<=x<=1/8. Denote L=ln(1/x).
The exact conditional energy moments imply
M1=E[R|cut]=a*x/(a+1),
Mlog=E[R(1+ln(1/R))|cut]
 =M1*(1+L+1/(a+1)).
Also
Mrem=E[R*|ln(x-R)| |cut]
 =M1*[L+psi(a+2)+EulerGamma].
For0<=a<=1 the last harmonic quantity is<=3/2.
At a=0 all these weighted quantities are zero; the possible
R=x endpoint is null for a>0. No expansion of F(a)x^a is made.


For a>0 the conditional energy density is a*r^(a-1)/x^a on(0,x).
Mlog follows by differentiating the exact power moment at p=1.
For Mrem substitute r=x*u and differentiate
Beta(a+1,b) at b=1. The digamma derivative is positive on the
positive axis, so psi(a+2)+EulerGamma<=psi(3)+EulerGamma=3/2
when0<a<=1. This is an analytic bound, not inferred from samples.

Independent quadrature tests both logarithmic integrals. Rare-index
tests retain a=10^(-800), and a separate resolution is10^(-10000).
The logarithm of the remaining energy is not silently replaced by
the logarithm of x: its extra beta-derivative term is positive.
