# Strict low-window moment and full-species comparison

Only the light stress cut is open for transfer T in [4,6] at
the actual mass, since 4M>6. Define

J_low=integral_4^6 rho_L(T)/(pi T^2) dT.

This is part of the form-factor slope, not the forward elastic
scattering cut calculated in S6.114.

## Controlled light-cut density on the window

Set K=(T-4)/4, y=x/(1-x), and u=M y(1+y). The parameter
density derived in cuts.md becomes

rho_L(T)=g/(16pi sqrt(T)) integral_0^K
 [y^2/((1+y)^3 M(1+2y))] du/sqrt(K-u).

For M>=1 and 0<=u<=K<=1/2, one has 0<=y<=1/2. Dividing
the square bracket by u^2/M^3 gives

1/[(1+y)^5(1+2y)].

The denominator lies between one and 243/16. Hence this ratio
is at most one and strictly above 1/16. The exact beta integral
is integral_0^K u^2/sqrt(K-u) du=16 K^(5/2)/15. Therefore,
strictly above threshold,

g K^(5/2)/(240pi sqrt(T) M^3)
 <rho_L(T)
 <g K^(5/2)/(15pi sqrt(T) M^3).

The lower estimate is used only on [5,6]. There K>=1/4,
sqrt(T)<3 and pi<4, giving rho_L>g/(92160 M^3).
On the full window K<=1/2, sqrt(T)>=2, pi>3 and sqrt(2)<3/2
give rho_L<g/(480 M^3). At T=4 the density is exactly zero.

The elementary inequalities 3<pi<4 are independently established
in the frozen ancestor chain by the arctangent integral and its
finite polynomial remainder; no decimal approximation to pi is used.

Since integral_5^6 T^-2 dT=1/30 and
integral_4^6 T^-2 dT=1/12, the strict positive moment obeys

g/(11059200 M^3) < J_low < g/(17280 M^3).

The endpoint strictness is supplied on open subintervals of
positive measure, not claimed for the zero threshold value.

## Comparison with the complete slope

The complete slope from S6.115 is

S=F_R'(0)=g/(96pi^2) integral_0^1 b(x)^2 dx,
b=x(1-x)/[xM+(1-x)^2].

For M>=1, the denominator is at most M, so
b>=x(1-x)/M. Its squared polynomial integral is 1/30.
With pi<4,

S>g/(46080 M^2),
0<J_low/S<8/(3M)<10^-196

at the actual rational mass. Also J_low<10^-603. Thus the
specified low-energy window has strictly positive but negligible
weight relative to the whole vertex slope.

The light and heavy full moments separately are

S_L=g/(96pi^2) integral (1-x)b^2 dx,
S_H=g/(96pi^2) integral x b^2 dx,
S_L+S_H=S.

Both are positive. The same lower bound on b gives the polynomial
integral 1/60 for each species, hence
S_L,S_H>g/(5760pi^2 M^2).
Using the already established upper bound S<g/(288pi^2 M^2)
gives S_L/S>1/20 and S_H/S>1/20.

This does not identify the remaining light moment with the
low window or identify M with a Regge scale. Both full moments
include the complete mathematical one-loop cut above their
respective thresholds. The full one-loop slope is not shown to
be dominated by the light threshold near invariant four.
