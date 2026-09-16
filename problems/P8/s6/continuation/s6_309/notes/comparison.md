# Uniform weighted variation and positive values of a named reference

Write P0(x) for the original elastic leading-soft reference.
For the marked state, the total emitted energy is R=omega<=x<=1/8.
S300/S301 establish, without expanding the detector power,

P_sigma(x)/P0(x) <= 1+4250/kappa.

The remaining-energy factor is then exact:

P_sigma(x-omega)/P0(x)
 = [P_sigma(x)/P0(x)]*(1-omega/x)^a_sigma
 <= 1+4250/kappa,

because a_sigma>=0 on the physical D4 state. This holds for arbitrarily
small positive x-omega; there is no lower remaining-energy cutoff and
no assumption that a_sigma*|ln(x-omega)| is small.

For a finite signed measure, the magnitude of its integral against a
nonnegative weight bounded by C is at most C times its total variation.
Applying that elementary inequality to the SAME S307 measure gives

|D[R](x)|/P0(x)
 < (1+4250/kappa)*min(1e32,2e14*x+2e37*x^2)/kappa.

The first factor is below 2 at the original kappa=1e800. This proves
the uniform constant bound 2e-768 and the uniformly vanishing bound
(4e14*x+4e37*x^2)/kappa. They hold at every nonforward angle without
a minimum transfer window, after the prescribed regulator removal.

Unlike adding the undressed remainder directly, the dressed comparison
does not need to divide its variation by the small power x^a_elastic:
that power is already in the state-correct convolution kernel.
No coupling expansion of that power is performed.

Since the relative correction is below one,

P0(x)+D[R](x) > P0(x)*(1-2e-768) > 0.

This is pointwise positivity for each hard configuration and threshold.
It does NOT prove that this approximation is monotone in threshold,
defines a positive measure on all radiation events, obeys unitarity, or
equals the full physical detector rate. The signed kernel and its
state-dependent prefactor require no such interpretation for this bound.

Regulator removal was proved at fixed positive x. The resulting reference
admits the uniform x->0 estimate above. This is not a proof that all
three dimensional-regulator, forward and detector limits commute for
the original interacting theory.
