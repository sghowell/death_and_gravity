# P8-S6.65: homogeneous finite local mass-response subtraction

Date: 2026-09-08. Original P8 remains OPEN.

Continue the restricted prepared-source, mass-only response of S6.64.
The physical metric, original clock, S6.55 state and S6.60 profiles
are unchanged. The source is homogeneous in space and smooth,
with compact time support strictly after u0=-1/2. This checkpoint
computes the local adiabatic subtraction and its specified finite
dimensional counterterm component, not the remaining nonlocal integral.

## Acceptance gates

1. Derive the moving-mode normalization, frequency and mass-readout
   variations independently from the physical D-dimensional Hamiltonian.
2. Linearize the actual WKB/Riccati recurrence through adiabatic order
   four, retaining the time derivatives of both mass-source profiles.
3. Perform the radial integrals and their first dimensional jets before
   setting D=3, retaining the D-1 transverse polarizations.
4. Recover the zero, two and four derivative homogeneous curved poles
   of S6.63 with the proper quadratic-action/response normalization.
5. Combine the radial finite part at mu=m with the frozen dimensional
   counterterm variation; recover a self-adjoint finite local operator.
6. Check the zero-derivative finite part against the already-frozen
   full vector potential, give exact compact coefficients and a
   continuous source-norm bound, pin sources and pass regressions.

The spatially constant source is a Fourier K=0 sector; no compact
spatial-support assertion is made for that individual sector.
The arbitrary-K bare kernel remains available at S6.64.

## Scope

This is the local term in the decomposition

    delta J_ren = integral(delta J_exact-delta J_ad,0..4)
                 + finite_local_component.

The remaining integral, its justified regulator limit, its quantitative
bound and its compatible-dimensional preparation estimates are still
required. The finite local component is evaluated in the stated
prescription with no additional finite Y^2 Wilson coefficient included;
it does not determine such matching coefficients from a UV theory.

The finite zero-derivative term agrees with the existing full potential
and is not to be added twice. All displayed source-quadratic derivative
terms vanish to first order on the selected clock. No profile or
background cancellation is changed.

The bound is a homogeneous C4-to-C0 source norm for the finite local
component only. It is not a contraction bound for the coupled equations,
nor a quantum stability, cone, interaction cutoff, V/G/B or P8 closure.
