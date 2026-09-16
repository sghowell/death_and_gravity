# One finite radiative residual with complete leading-soft dressing

## Original source and named operation

Retain mu=nu=1, kappa=1e800, 25/4<=s<=16, every physical nonforward
hard angle and 0<x<=1/8. The positive complete Born amplitude is
A_B=A_m+A_G. S295/S304 own all 26 matter and 21 Einstein radiation
trees; S307 owns their finite signed measure

dR(sigma)=[rho sum_pol |M5/A_B|^2-sum_pol |S0|^2] dPhi1.

Here sigma is the exact one-graviton radiative state, of marked energy
omega, with the specified massive recoil. Its cumulative total variation
obeys V_R(x)<min(B,L*x+Q*x^2)/kappa, with
B=1e32, L=2e14 and Q=2e37. No graph or interference is discarded.

Define only the following SINGLE-RESIDUAL reference:

D[R](x)=integral_(0<omega<x) P_sigma(x-omega) dR(sigma),
P_sigma(y)=exp(Delta_sigma)*F(a_sigma)*y^a_sigma,
F(a)=exp(-gamma_E*a)/Gamma(1+a).

The full conserved massive-plus-marked-null current fixes a_sigma and
Delta_sigma in the unchanged analytic convention. The marked hard state
is held fixed while any number of additional LEADING-soft emissions is
summed. This operation is not asserted to be an exact decomposition of
the complete interacting multi-real rate.

## Regulator limit and main result

At soft regulator e>0 the additional-N term, including its Bose factor,
is [a_e Gamma(2e)y^(2e)]^N/[N! Gamma(1+2eN)].
Its virtual factor is exp[-a_sigma/(2e)], not the elastic virtual pole.
The fractional-D angular contraction need not be positive; a uniform
absolute series bound, rather than an unsupported probability law,
justifies integrating the limit against the signed residual.

At every fixed positive detector threshold this sum converges to the
displayed P_sigma, including the a_sigma=0 boundary. The total-energy cut
requires y=x-omega. The Bose relabeling argument supplies exactly one
marked residual and no extra factorial or independent-energy cut.

For the elastic reference P0(x), the complete dressed correction obeys

|D[R](x)|/P0(x)
 < (1+4250/kappa)*min(B,L*x+Q*x^2)/kappa
 < 2e-768,

and also

|D[R](x)|/P0(x) < (4e14*x+4e37*x^2)/kappa -> 0.

Both estimates are uniform in every nonforward hard angle. The powers of
the detector and remaining energies are never expanded. P0+D[R] is
therefore positive pointwise on this domain.

## Strict boundary

Pointwise positivity of this named reference does not prove monotonicity
in x, a positive detector measure, unitarity or equality with the complete
physical rate. Multiple finite residuals, simultaneous two-soft contact
terms, correlated exact recoil, finite radiative hard loops and physical
matching require separate estimates. Original V/G/B/P8 remain OPEN.
The dimensional regulator is removed first at fixed positive x; no
arbitrary joint dimensional/forward/resolution limit is asserted.
