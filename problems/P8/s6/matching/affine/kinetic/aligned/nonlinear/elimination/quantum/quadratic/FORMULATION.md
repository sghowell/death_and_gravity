# P8-S6.61: curved quadratic retained-vector mass-insertion pole

Date: 2026-09-08. Original P8 remains OPEN.

The S6.60 selected-state action, vector state and clock background are
unchanged. This checkpoint derives the missing four-dimensional local
quadratic mass-insertion UV residue. It does not evaluate the finite
retarded second variation or replace S6.53's finite prescription.

## Objects and conventions

Work with the Euclidean generalized Proca operator on a positive
background metric and a constant reference mass m>0,

    M^{mu nu}=m^2 g^{mu nu}+Y^{mu nu}.

The displayed quadratic density Q obeys

    Gamma_div^(Y^2)=(32 pi^2 epsilon)^(-1) integral sqrt(g) Q,
    epsilon=(4-d)/2.

The Fourier coefficients are bilinear in opposite spatial momenta.
All integration-by-parts comparisons concern compactly supported
variations, not deletion of asymptotic boundary fluxes of the bounce.
The Lorentzian expression requires the same Wick/sign conversion as
the frozen dimensional matching calculation.

## Acceptance gates

1. Derive ordered tensor derivatives on a Euclidean FLRW frame and
   independently derive coordinate curvatures of the auxiliary metric.
2. Obtain orders zero and two from the covariant mass-insertion basis;
   check order two against the independent covariant-symbol calculation.
3. Derive all order-four mass dependence from the scalar determinant
   on g_tilde=sqrt(det(1+Y/m^2)) (1+Y/m^2)^(-1) g. Expand the full
   scalar heat coefficient, without using a four-dimensional
   Gauss-Bonnet reduction to define the operator basis.
4. Replay all three flat Feynman-bubble coefficients; check the
   curved scalar-mass conformal action and constant anisotropic
   mass scaling.
5. Insert the actual first lapse derivatives of the retained mass
   tensor, including their time derivatives, and derive an exact
   compact-support self-adjoint local action on the original clock.
6. Retain the failing literal fourth-order reference transcription
   as an explicit withheld control, not an accepted pole.
7. Pin sources, rebuild S6.60, reject unsupported inputs and pass
   independent ordinary and full P8 regressions.

## Matching and scope

S6.53 already includes the full flat potential, including its Y^2
term. That term is checked, not added a second time. The new
two- and four-derivative Y^2 operators have zero value and first
variation at Y=0. They therefore supply missing quadratic pole
data without changing the selected clock tadpole or S6.60 profiles.

This checkpoint supplies the four-dimensional residue and an
explicit covariant operator definition. A chosen dimensionally
continued counterterm, its finite evanescent contributions, the
full mixed metric/mass second variation, compatible varied-state
preparation and the integrated renormalized retarded kernel still
require calculation. No pole residue is treated as a finite error
bound, exact higher-derivative degree-of-freedom count, quantum
cone, cutoff or V/G/B result. Original P8 is not closed.
