# Actual map allowance and full analytic target comparison

At kappa=10^800 the exact finite gate majorant is
approximately 4.2172416e-6392 and is strictly below
the declared rational cap 10^-6390.
The unchanged cubic coefficient norm is below
10^-403, so the map difference is strictly below
the declared delta=10^-6793.

The actual M is far above 4356. The exact
linearized action allowance is approximately
3e-6791 and strictly below 10^-6790. This
includes all free, local quartic and full
heavy-resolvent terms.

S6.162 already bounds the old mapped polynomial
stationary action against the FULL finite-kappa
analytic target on this same class. Add the
new map allowance by the triangle inequality:

    |S_polynomial[Fhat(Psi),H_stationary]
      -S_full_analytic_target[Psi]|
    <=(E_previous+E_map_change) ||Psi||_2^2
    <10^-800 ||Psi||_2^2.

The final inequality is evaluated on exact
rationals, not by rounded decimal comparison.
The SAT8 mass modification vanishes identically
when the classical fermions are zero, so it
does not alter this stationary-action comparison.

This proves neither quantum matching of the
analytic target nor an estimate of physical
omitted loop orders. All boundaries remain
Schwartz and the full inverse is used only on
the explicitly stated Fourier class.
