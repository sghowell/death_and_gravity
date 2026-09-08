# Exact massive integral bounds and the noncancelled pole

Put gamma=1/(64*pi^2*L^2) and define, for real s>0,

Kfr(s)=E0+s^2 E2+gamma*(s^4 Breg(s)+C0+s^2 C2+Ctad),
Breg(s)=Freg-(1/4) integral_0^1 y^2 Mreg(y^2)/(d-y^2) dy,
d=1+4*m^2/s^2.

The positive pair matrix is Mreg(z)=2 bT bT^T+bL bL^T,
with
bT=((109/81)*(1-z),2*(1-z)),
bL=((109+117*z)/81,2*(1+z)).
These follow from the actual physical Hamiltonian vertices.

All pair entries are nonnegative for 0<=z<=1. Their endpoint
upper bounds give, entrywise,

Mreg <= [[24946/2187,1340/81],[1340/81,24]]
     <= [[12,17],[17,24]].

Since Freg_NN=-7357/6561<-1 and the integral is positive
semidefinite, Breg_NN<=Freg_NN<-1 for every positive s.

## A continuous, not sampled, band bound

For 0<=y<=1 and d>1,

y^2/(d-y^2) <= 1/(d-y),

because d-y^2-(d-y)=y*(1-y)>=0. Hence its integral is at most
log(d/(d-1))=log(1+s^2/(4*m^2)).

For m=1000 and s<=10^13 the logarithm is below 100. This
uses the exact positive exponential series: its single
term 100^18/18! exceeds
1+10^26/(4*10^6)=25000000000000000001.
Together with abs(Freg) <= [[2,3],[3,4]], this proves

abs(Breg_ij) <= [[302,428],[428,604]]_ij

throughout the whole closed band [10^12,10^13].

The elementary bounds 9<pi^2<10 need no floating-point
evaluation. The first eight alternating terms of
4*integral_0^1 1/(1+x^2) dx give 135904/45045>3, with
positive remainder integrand 4*x^16/(1+x^2).
The positive integral of x^4*(1-x)^4/(1+x^2) on [0,1]
equals 22/7-pi, and (22/7)^2<10.

All lower-order terms therefore satisfy, entrywise and
uniformly on the band,

abs(gamma*(C0+s^2 C2+Ctad))
 <= (2*10^15+10^26*10^9)/(576*L^2)
 < e, e=10^-12.

The exact (smaller) rational upper bound is in the report.

## Determinant signs, allowing nonsymmetric lower rows

Let c=749377/250000, g-=s^4/(640*L^2),
g+=s^4/(576*L^2). Write the four actual entries of Kfr
as A,B,C,D in row-major order. The following bounds hold:

A >= a-=c-302*g+-e,
A <= a+=c-g-+e,
D >= d-=6*s^2-9/100-604*g+-e,
abs(B),abs(C) <= o=3/200+428*g++e.

At s=10^12, exact rational arithmetic gives
a->2, d->5*10^24, o<1. Thus
det Kfr=AD-BC >= a-*d--o^2 > 0.

At s=10^13 it gives
a+<-12, d->5*10^26, o<10000. Since a+<0 and D>=d->0,
AD<=a+*d-, and therefore
det Kfr <= a+*d-+o^2 < 0.

These inequalities do not require B=C. Both endpoint
determinant comparisons are checked as exact rationals.
Moreover the band-wide lower bound obtained by replacing
s^2 by its lower endpoint and g+ by its upper endpoint
keeps D strictly positive everywhere in the band.

The exact massive bubble is analytic for Re s>0 and real
continuous for positive s, by its positive-mass integral
or S6.72's first-sheet analysis. The local terms are
polynomials and the tadpole values are constant. Thus the
intermediate value theorem supplies at least one real
s0 in the open band where det Kfr(s0)=0.

The determinant is not identically zero, since both
endpoints have strict signs. Its analytic zero is isolated
and has finite multiplicity. The NN entry of adj(Kfr) is
D(s0)>0, so (Kfr^-1)_NN=D/det Kfr has a nonremovable pole.
No assumption of a simple zero is needed.

The tree-only determinant is
6*c*s^2-9*c/100-9/40000, strictly positive throughout this
band. S6.72's isolated Breg is also invertible there.
It is their full coupling that has the certified pole.

## What this says about stability

A bounded causal convolution transfer on the entire
half-line has a bounded analytic Laplace transform in
Re s>0 (by its instantaneous-plus-integrable kernel
representation). It cannot equal Kfr^-1 near this pole.
The frozen continuum transfer consequently has no such
bounded stable inverse. A finite-interval causal inverse,
whose norm may be enormous, is perfectly compatible with
this obstruction; S6.73 did not assert half-line stability.

The result is about this explicitly defined frozen model.
There is no stationary Laplace transform of the actual
time-dependent operator to which the intermediate value
argument can simply be reapplied. In particular, this
proof neither bounds the omitted curved/state kernel nor
establishes an actual bounce growth rate.
