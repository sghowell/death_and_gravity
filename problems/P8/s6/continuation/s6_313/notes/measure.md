# Complete signed two-real measure

Define the density by |F2|^2-|B2|^2, not by |R2|^2 alone.
The identity
 |F2|^2-|B2|^2=2Re(B2*conjugate(R2))+|R2|^2
retains every interference. Each unsubtracted integral diverges at
soft endpoints; the difference is defined with a common lower cutoff.

Put a=W*z,b=W*(1-z), I=W*h(z),
h(z)=-z ln z-(1-z)ln(1-z), da db=W dW dz.
The exact sharing integrals are
 int h(z)/[z(1-z)] dz=pi^2/3,
 int h(z)^2/[z(1-z)] dz=4*zeta(3)-pi^2/3<2.
For the first, expand -ln(1-z) and integrate its positive series.
For the second, the two squared-log terms each contribute
2*(zeta(3)-1). The cross term uses
sum_(n>=1)1/[n(n+1)^2]=2-zeta(2), obtained from
1/[n(n+1)^2]=1/n-1/(n+1)-1/(n+1)^2.
Monotone convergence justifies each positive termwise integral.
zeta(3)<5/4 follows from the decreasing-series integral bound,
and zeta(2)>3/2 from a finite positive partial sum; hence the strict
upper bound2. Independent70-digit quadratures agree within10^-60.

Using |R2|<C*I/(ab), |B2|<D/(ab), the radial integrals obey
 int ab*2|B2 R2| da db<(2*pi^2/3)*D*C*x<7DCx,
 int ab*|R2|^2 da db
 <C^2*x^2/2*[4*zeta(3)-pi^2/3]<C^2*x^2.

Each graviton phase measure is
a da dOmega/[2(2pi)^3]. The two-identical-graviton factor1/2!
gives1/[8(2pi)^6] before angular integration. Summing the four
unit physical polarization pairs and bounding both solid angles
by4pi yields1/(8pi^4)<1. The rho2 factor was already included in F2,
and the same phase convention was built into both face terms.
Therefore the total variation of the whole signed measure is
 <7DCx+C^2x^2<2*10^-725*x+10^-652*x^2.
This integrable majorant yields the cutoff-independent finite signed
measure by dominated convergence, including all angular endpoints
in the bounded almost-everywhere sense. No choice of collinear
polarization coordinates changes the absolute bound.

The frozen unexpanded elastic leading-soft reference obeys
P0>x^alpha/2 with0<=alpha<1/2. For0<x<=1/8, division gives
TV/P0<4*10^-725*sqrt(x)+2*10^-652*x^(3/2), uniformly tending to zero.
This is a comparison of size only. It does not identify |B2|^2 with
previous virtual or real subtraction terms, nor prove positivity,
normalization or monotonicity of a complete detector distribution.
