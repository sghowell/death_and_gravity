# Flat full-tensor absorptive benchmark

## Fixed object

Retain the S6.176 canonical ordinary Proca sector, mass m=1000,
kappa=10^800, and the full parent normalization. Compute its connected
stress two-point function in a flat reference vacuum. This reference
calculation does not replace the actual CD state, choose a new finite
prescription or constitute a complete parent vacuum scattering amplitude.

The observable is the external-metric retarded current
i theta(t-s) <[T_ab(t),T_cd(s)]>/4. Retain the full symmetric stress,
three physical polarizations, constrained temporal readout, all nine
creation pairs, angular tensor average, Wick factor and Lorentz measure.

## Statement

For s above 4m^2 and beta=sqrt(1-4m^2/s), the conserved spin densities are

    rho2(s) = beta (13s^2+56m^2 s+48m^4)/(3840 pi^2)
    rho0(s) = beta (s^2-4m^2 s+12m^4)/(384 pi^2).

Both are zero through threshold and positive above it. The complete
timelike tensor is rho2 P2+rho0 P0; P0 has rank 1 and P2 rank 5 in
the three-dimensional transverse space.

The dispersive PART

    D_i(z) = z^3 integral_[4m^2,infinity) rho_i(s)/(s^3(s-z)) ds

is analytic off the physical cut. For |z|<=2m^2,

    |D_i| <= 2 |z|^3 I_i
    |D_i-z^3 I_i| <= 2 |z|^4 J_i,

where

    I2=3/(3584 pi^2 m^2), I0=3/(8960 pi^2 m^2),
    J2=31/(322560 pi^2 m^4), J0=1/(32256 pi^2 m^4).

Explicit inverse-spectral-limit tail bounds hold. No physical polynomial
or counterterm is inferred from this absorptive information.

The four-dimensional vector-minus-scalar heat weight is

    13 C^2/120 - 7 Euler/40 + R^2/72

modulo a separately tracked total derivative. The flat C^2 and R^2
Hessians supply the respective factors s^2 P2 and 6s^2 P0. These agree
with the ultraviolet densities, including the longitudinal contribution.
This is not dimensional continuation of the finite local action.

Both external-metric changes of normalization multiply the response
by 4/kappa. The spin0 tensor projector is not a canonically reduced
propagating scalar or a full constraint-reduced mixed norm. The spectral
density still grows as s^2, so no uniform all-energy smallness follows.

## Proof and boundary

The notes derive physical amplitudes, complete projection and current
normalization, dispersion and limits, curvature comparison, scope, and
independent validation. Exact identities certify the finite algebra;
the analytic and continuous arguments are written proofs, not FORMALIZED.

The S6.198 full reference contact and first five equal-time kernels,
their fixed covariant spatial matching, response inverse, nonlinear
background and stability remain unresolved. So do the other original
V/G/B obligations. Scoped P8(a) and A.20-A.23 are unchanged.
