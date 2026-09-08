# P8: the local nonlinear auxiliary-constraint gap is closed

Original P8 remains open. The scoped photon objective, linear CD/M1
classification, physical matter frame and adopted V/G/B contract are
unchanged. This follows the [on-shell linear-light checkpoint](assessment-2026-09-07-p8-onshell-light-family.md).

## New result

[S6.45](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/FORMULATION.md)
derives the complete scalar-clock ADM and Hamiltonian constraints of
the unchanged source-aligned S6.42 action. Keeping both the spatial
curvature boundary and the linear lapse-velocity primitive removes
all lapse time and spatial derivatives. The original scalar potential,
full retained source/mass, Maxwell field and free chi are retained.

The trace/temporal Legendre calculation gives effective gamma_t>9/10
on the original closed clock tube. All five shear, three vector and
free-matter velocity pivots are included. Performing the spatial Maxwell
Legendre transform before using the normal temporal variable shows
that its canonical equation is algebraic; the one-form spatial
diffeomorphism generator includes its temporal Gauss term.

On the actual rolling solution, the complete lapse/normal-vector
auxiliary Jacobian divided by sqrt(h_hat) is diag(-2J,-1).
The original J>1/40 throughout [-1/2,1/2], including the bounce,
without a Theta divisor. Smoothness and compactness give a nonempty
local canonical/spatial-jet neighborhood with an invertible auxiliary
Jacobian. This is an existence theorem for that neighborhood, not a
numerical field/jet radius or a theorem on the entire X tube.

The four auxiliary constraints have a local block inverse even when
their secondary-secondary brackets contain spatial differential
operators. Their consistency fixes the multipliers. Together with
the six first-class spatial-diffeomorphism/shift constraints, the count
is seven physical modes including free chi: two tensor, one clock
scalar, three vector and one chi. The all64 complement and projective
gauge directions were already correctly removed by S6.41.

## What this closes, and what remains

This closes the local nonlinear secondary-constraint/count gap for
this unchanged candidate. It does not establish nonlinear energy,
stability or causality on other backgrounds, global well-posedness,
the X=0 vacuum count, full nonlinear light/heavy evolution, higher-order
or quantum errors, an interacting cutoff, or V/G/B UV admissibility.

Work continues on the exact fixed-light vector elimination and the
full nonlinear source remainder. The retained vector is Gaussian at
fixed light fields, but a formal inverse is not a numerical remainder
bound; retarded response is not a one-copy variational action. The
quantum determinant is not deleted by classical source alignment.
No new user decision or authorization is currently required.

## Verification

The report pins 16 sources and fully rebuilds S6.44 and its ancestry.
It checks 23 named exact identities comprising 47 scalar entries,
11 rank/continuous/interface proof checks and 7 rejected inputs.
The independent scientific suite passes **32 tests in 3.91 seconds**;
the ordinary suite passes **53 tests in 184.96 seconds** with no
broad GCD adapter. The separate seeded read-only CLI passes.

The complete P8 regression through S6.45 passes **4497 tests in
859.57 seconds**, with no frozen checkpoint excluded. The unchanged
exact GCD runner passes its 128 original tuple controls and records
5706 domain fallbacks and 6509 exact descents, using the existing
per-test seeded recipe and no faulthandler timers.

Before freezing, the controls exposed an import-root typo, an inexact
test literal and an installed SymPy block-collapse simplification
defect. The latter was isolated and avoided by materializing the exact
block product before identity subtraction; a regression control remains.
No mathematical criterion or frozen ancestor was weakened or changed.
Verification is exact symbolic checking with written continuous and
Dirac-closure proofs, not proof-assistant formalization or external review.
