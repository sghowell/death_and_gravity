# Relative comparison without expanding the resolution power

The named elastic leading-soft total-energy reference is

P(x)=exp(Delta)*F(alpha)*x^alpha,
F(alpha)=exp(-gamma_E*alpha)/Gamma(1+alpha).

S299 proves the complete all-N leading-soft sum with the total emitted
energy cut, rather than a product of individual energy cuts.
On the original compact window its soft index is in fact below
4/(5kappa); the looser0<=alpha<6/kappa suffices here.
S296 gives|Delta|<112/kappa.
At originalkappa=1e800 these implyalpha<1/2 and|Delta|<1/4.
The Weierstrass product givesF>=1-alpha^2>3/4.
Alsoexp(Delta)>=1-|Delta|>3/4. HenceP(x)>x^alpha/2.

Letx0=1/1000. Forx<=x0 bothx^(1-alpha) andx^(2-alpha) increase.
Thus the new absolute bound divided byx^alpha is at most its value
atx0. The numerator there is

2e14*x0+2e37*x0^2
=20000000000000000000200000000000<1e32.

Forx>=x0 use S304's independent fixed1e32 bound and
x^-alpha<=x0^-alpha. Both cases give

TV(x)/x^alpha<1e32*1000^alpha/kappa.

The exact finite positive Taylor sum through degree8 provesexp(9)>1000,
so ln1000<9. Sincealpha<6/kappa and54/kappa<1/2,
1000^alpha<exp(1/2)<2; the last inequality follows by bounding the
exponential series by the geometric series, strictly from degree2 on.
Together withP(x)>x^alpha/2,

TV(x)/P(x)<4e32/kappa=4e-768

for every0<x<=1/8 and every nonforward hard angle. This remains valid
for exponentially fine resolution wherealpha*lnx is not small.

To make the uniform zero-threshold conclusion explicit, usealpha<=1/2
and0<x<1 directly in the new bound:

TV(x)/P(x)<[4e14*sqrt(x)+4e37*x^(3/2)]/kappa ->0.

This is a comparison theorem, not an all-N nonleading error theorem.
Adding the single-real remainder toP(x) defines a useful positive
known approximation (because its relative correction is tiny), but
does not prove equality or a remainder bound for the complete detector
observable. No such promotion is made in the certificate.
