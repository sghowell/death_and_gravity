# Continuous coefficient error budget

S6.176 gives |r^(j)|,|p^(j)|<epsilon=1e-770 for j=0,...,5
on I=[-1/2,1/2]. Also 0<p_v<1e-789. These are absolute
normalizations by fixed kappa, not by the bounce energy,
which vanishes.

## Full original X strip

For every real x in (-1/4096,6/5), take a complex disc
|z|<=d=1/(64N)=1/65536 about x. T is holomorphic and
|T|<=2 on every such disc, as follows.

If x<=1/4, the ratio |(x+z)/(1-x-z)| is at most
(1/4+d)/(3/4-d)<1/2. Write T=q^N/(1+q^N).
If x>=3/4, the reciprocal ratio obeys the same bound;
write T=1/(1+q^N). These estimates exclude denominator
zeros and give |T|<2.

For 1/4<=x<=3/4, both real parts x+z and1-x-z
are at least1/4-d>=1/8. Since arctan(y)<=y for y>=0,
the sum of their absolute arguments is at most16d.
The ratio q=(1-x-z)/(x+z) therefore has
|N arg q|<=16Nd=1/4<pi/2. Its Nth power has positive
real part, so |1+q^N|>=1 and |T|<=1. This covers the
transition where a naive real denominator estimate would
be far too large.

Cauchy's formula gives the exact outward envelope

    |T^(k)(x)|<=2 k! (64N)^k, k=0,...,5.

On this strip |X-1|<2. The bracket V=-p+p_v+B(X-1)
has |partial_u^j V|<4epsilon and
|partial_u^j B|<epsilon for j<=5.
For X order0, |partial_u^j DeltaF|<5epsilon.
For1<=k<=5, the product rule yields

    |partial_u^j partial_X^k DeltaF|
       <[4*(2 k! (64N)^k)
          +k*(2 (k-1)! (64N)^(k-1))]epsilon.

Every one of the36 mixed derivatives is below1e-742.
The rational Cauchy constants and all36 comparisons are
checked exactly. This budget includes the null-gradient
locus and the connected X interval between vacuum and clock.

## Sharper clock-domain budget

For3/4<=x<=6/5 choose the larger complex radius1/32.
On its disc |(1-x-z)/(x+z)|<=(9/32)/(23/32)=9/23<2/5.
Thus |T-1|<2(2/5)^1024, and Cauchy gives

    |(T-1)^(k)|<2 k!32^k(2/5)^1024<1e-390

for k=0,...,5. For real x, 0<=T<=1 as well.
Here |X-1|<=1/4 and |partial_u^j V|<2epsilon.
The36 mixed derivatives of DeltaF are each below2epsilon.
For X orders2 through5 they are, more sharply, below
10epsilon*1e-390=1e-1159. These estimates do not require
real analyticity of the reference stress in u.

## The actual local lapse pivot

Write delta=1/[2(1+u^2)^3]. The literal full chart has
R(N)=1+2delta(N^-2-1), volume R^-3/4 and X=N^-2.
The addition to the classical Hamiltonian lapse pivot is

    DeltaJ=(21delta^2-3delta)A/2+(1-6delta)B.

This follows by twice differentiating N R^-3/4 DeltaF
at N=1; it is not a dropped-boundary clock substitution.
For32/125<=delta<=1/2 its A coefficient is between0 and
15/8, and the B coefficient has absolute value at most2.
Thus |DeltaJ|<4epsilon.

The actual old J has a positive-coefficient numerator with
constant1215 and denominator800(1+u^2)^18. Therefore on I

    J_new>=1215/[800(5/4)^18]-4epsilon>1/100.

This preserves this specific classical algebraic constraint
pivot, not the full quantum response spectrum or stability.
