# Signed dimensional sum and its uniform dominating function

Let a=a_sigma>=0 be the physical D4 soft index and let

a_e=p(e) K_e/(4*pi^2*kappa), 0<e<=1/8.

The unchanged angular convention gives
a_e=a+2e*Delta_sigma+O(e^2). The full-D angular contraction on the fixed
three-ball is real but is not assumed positive in the interior.
Consequently a_e is not assumed nonnegative uniformly over all states.

## A uniform bound without a finite-e probability assumption

S301 gives |K0|<265, |K1|<130000, |p'(0)|<4 and

|p(e)K_e-K0-e(K1+p'(0)K0)| < 2000000*e^2.

Using pi>3 and e<=1/8 therefore yields

|a_e-a|/(2e)
 < [130000+4*265+2000000/8]/(72*kappa)
 = 10585/(2*kappa) < 5500/kappa.

At remaining energy 0<y<=1 define the virtual-normalized analytic sum

S_e(y)=exp[-a/(2e)] sum_(N>=0)
 [a_e Gamma(2e)y^(2e)]^N/[N! Gamma(1+2eN)].

For every z>=0, Gamma(1+z)>=integral_1^infinity exp(-t)dt=exp(-1);
hence its reciprocal is at most exp(1). Log-convexity between the
endpoints Gamma(1)=Gamma(2)=1 implies Gamma(1+2e)<=1. With
Gamma(2e)=Gamma(1+2e)/(2e) and y^(2e)<=1, the ABSOLUTE sum is at most

exp(1)*exp[(|a_e|-a)/(2e)]
 <= exp(1)*exp[|a_e-a|/(2e)]
 < exp(1)*exp(5500/kappa) < 6

at the original kappa. The last step uses exp(1)<3 and
5500/kappa<1/2, so exp(5500/kappa)<2.
All factorial sums are absolutely convergent. This bound is uniform in
the radiative state, y, multiplicity and the admitted regulator.
It needs no positivity of a_e and no reinterpretation of fractional
dimensions as physical helicities.

## The two distinct physical-index cases

At any fixed state with a>0, a_e is eventually positive. Set
lambda=a_e Gamma(2e)y^(2e). The positive series can then be written

S_e=exp[lambda-a/(2e)] E[1/Gamma(1+2eN)],
N~Poisson(lambda).

The auxiliary scaled variable obeys 2eN->a in L2, while

lambda-a/(2e) -> Delta_sigma-gamma_E*a+a*ln(y).

The reciprocal Gamma function is bounded and continuous on the positive
axis, so the expectation tends to 1/Gamma(1+a). This is a limit of the
whole sum, not a limit taken at each fixed N.

If a=0, lambda instead tends to Delta_sigma, which may have either
sign in the abstract analytic continuation. For e in a sufficiently
small interval, |lambda| is bounded. The terms are dominated by
exp(1)*C^N/N!, and 1/Gamma(1+2eN)->1 at every fixed N. The limit of
the entire signed series is exp(Delta_sigma), including negative Delta.
This agrees with the same formula at a=0.

Thus at fixed y>0 both cases give
P_sigma(y)=exp(Delta_sigma-gamma_E*a)*y^a/Gamma(1+a).
The common bound by 6 permits dominated convergence against any finite
signed marked measure of the present type. In the physical residual,
omega=x is a zero-measure endpoint, so no pointwise y=0 theorem is needed.
No arbitrary simultaneous e*|ln(y)| scaling is claimed.
