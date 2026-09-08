# P8-S6.62: quadratic counterterm dimensional continuation

Date: 2026-09-08. Original P8 remains OPEN.

This continuation preserves the S6.60 selected action and clock
profiles and all lower-order finite matching choices frozen in S6.53.
The local expressions retain the Euclidean conventions of S6.61;
their Lorentzian Wick/sign conversion is a separate required step.
Only the previously unspecified two- and four-derivative quadratic
mass-deviation counterterms are considered. The full flat potential
is not added again.

## Literal continuation

Use D=3-2epsilon spatial dimensions. Continue every contraction of
the S6.61 two-derivative quadratic basis covariantly in D+1 dimensions,
with its four-dimensional numerical coefficients held fixed. The
temporal and common spatial eigenvalue deviations remain scalar
coefficient functions independent of D.

For the four-derivative term keep the literal covariant definition

    g_tilde,D(lambda)=sqrt(det(1+lambda Y/m^2))
                     (1+lambda Y/m^2)^(-1) g_D,
    Q4,D=-[lambda^2]sqrt(g_tilde,D)/sqrt(g_D) a4sc(g_tilde,D),
    a4sc=R^2/72-Ricci^2/180+Riemann^2/180.

The determinant and tensor ranks are evaluated in D+1 dimensions.
The determinant exponent 1/2 and the curvature coefficients are
fixed. This is a specified counterterm continuation, not a claim
that the auxiliary metric's scalar-action interpretation is
unchanged away from four dimensions. Other evanescent choices
would be different finite prescriptions.

The counterterm sign is -Q_D/(32 pi^2 epsilon). Define its
normalized lapse variation E_D=a^(-D) delta(integral a^D Q_D)/delta n.
The pole and evanescent parts are

    -E_3/(32 pi^2 epsilon)+2(partial_D E_D)_3/(32 pi^2).

Differentiate and integrate by parts with the actual D-dependent
measure before taking the dimensional limit.

## Acceptance gates

1. Derive generic-D coordinate curvature from a flat warped fiber
   and check every invariant against the independent S6.61 geometry.
2. Evaluate all two-derivative Einstein sums by exact transverse
   index partitions, not numerical dimensional interpolation.
3. Recover every S6.61 pole sector at D=3 and derive its first D jet.
4. Insert the actual time-dependent lapse mass profiles, derive the
   D-dependent compact local action and compare its variation with
   direct variation of the unreduced density through first order in D-3.
5. Display exact first-dimensional-jet local coefficients and retain
   the measure-variation terms in the normalized response operator.
6. Prove a nonzero finite Gauss-Bonnet boundary control, pin sources,
   fully rebuild S6.61, reject unsupported inputs and pass regressions.

## Scope

Each added quadratic operator has zero value and first variation
on the clock in every continued dimension. No existing selected
stress profile is changed. The result specifies and computes an
evanescent part of the new counterterm; it is not the finite
renormalized response by itself.

The dimensionally continued bare kernel, compatible varied-state
data, its contact/retarded terms and common normalization must
be combined with the counterterm before a finite limit. No finite
Wilson-coefficient boundary values or integrated feedback bound
are supplied here. Coupled quantum stability, cones, interactions,
cutoff and V/G/B remain open.
