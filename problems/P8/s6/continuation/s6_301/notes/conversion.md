# Finite soft reference conversion and its remaining hard boundary

On the SAME recoiled radiative state, the universal real pole is
K0/(8*pi^2*kappa*e); the virtual rate pole has the opposite sign.
Using an elastic virtual pole with a radiative real state instead leaves
a genuine (K0_state-K0_elastic)/e divergence. Its small numerator does
not make it finite at fixed positive radiation energy.

With the inherited analytic soft convention and fixed resolution x,
the universal paired soft factor is

x^(2e)*[p(e)*K_e-K0]/(8*pi^2*kappa*e).

Its finite conversion is

Delta_sigma=[K1+(EulerGamma-2-lnpi)*K0]/(8*pi^2*kappa).

The coefficient bounds, |EulerGamma-2-lnpi|<4 and pi>3 imply
|Delta_sigma|<2000/kappa. The explicit second-order remainder and
1-x^(2e)<=2e*|lnx| show

|paired(e)-Delta_sigma|
 <= e*[28000+4000|lnx|]/kappa

for e in [0,1/8], x in (0,1/8], nu=1. Remove e at fixed x.
This is not a uniform limit for arbitrary e*|lnx| scaling.

## Known detector-reference factors only

Let F(a)=exp(-EulerGamma*a)/Gamma(1+a). Compare the known functions

exp(Delta_sigma)*F(a_sigma)*x^a_sigma

and the corresponding elastic expression, only after correct same-state
pole pairing. S300 already supplies |delta log(x^a F(a))|<12/kappa+
60/kappa^2 when R<=x<=1/8. The present bound and S296 give
|Delta_sigma-Delta_elastic|<2112/kappa. At original kappa=10^800,
the complete log ratio is<2125/kappa<1/2, hence the multiplicative
ratio differs from one by<4250/kappa<5*10^-797.

Separately, a_sigma<6/kappa and the Gamma-product bound give
0<=1-F(a_sigma)<=a_sigma^2. At original parameters |Delta_sigma|<1/2,
so |exp(Delta_sigma)-1-Delta_sigma|<=Delta_sigma^2 and
exp(Delta_sigma)<2. Therefore

|exp(Delta_sigma)*F(a_sigma)-(1+Delta_sigma)|
 < 5000000/kappa^2 = 5*10^-1594.

The detector power x^a is NOT expanded. These are comparisons of
known leading-soft reference factors, not a bound on the full rate
or the all-multiplicity non-leading remainder.

## Scheme and evanescent hard terms

The angular/projector/phase prescription is fixed here. If a hard
tree or loop is A0+e*A1+..., its A1 times an infrared pole can
contribute a finite hard term. Neither K1 nor Delta_sigma sets A1 to
zero. A full physical comparison must use the hard amplitude and
external-state continuation in the same scheme. This calculation
does not promote one universal factor to a scheme-independent hard
amplitude, complete radiative unitarity or a Regge theorem.
