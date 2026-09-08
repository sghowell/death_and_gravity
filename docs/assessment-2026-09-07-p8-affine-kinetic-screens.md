# P8: exact matching survives; these kinetic promotions do not complete it

Original P8 remains open. The completed scoped photon objective and
linear CD/M1 classification are unchanged. This follows the
[exact auxiliary lift](assessment-2026-09-07-p8-affine-cd-lift.md), not a
revision of its source-hashed action or proof.

## What the new calculation decides

[S6.38](../problems/P8/s6/matching/affine/kinetic/FORMULATION.md)
tests three explicit kinetic additions to the affine action.

The homothetic-curvature square gives exactly the original CD/M1 theory
plus a gravitationally coupled Maxwell field. The connection quotient
that carries the matching remains auxiliary, so this is a spectator
promotion and supplies no new UV control.

The simple matrix-trace square of trace-free curvature has opposite
time-kinetic signs in its independent spatial spin-3 and axial
connection sectors. Neither is removed by a projective, temporal or
metric constraint. No overall coefficient sign makes this named action
healthy. The full trace-free projector is retained and tested before
restricting to those traceless physical sectors.

The selective projectively invariant quotient-vector curl is more
instructive. Its remaining 56 connection variables can be eliminated
exactly, with a continuous four-vector inverse bound on the original
closed tube. All 64 original connection source components vanish on
the actual rolling background. The extra term therefore preserves the
same complete physical bounce and free-matter backreaction exactly.

Nevertheless, its coupled perturbations fail the kinetic-health test.
The full rolling vector source is derived before discarding any
metric, shift or coefficient variation. A true exact-gradient change
removes the apparent new lapse-time derivative; thus the failure is
not a mistaken claim that a raw lapse-velocity square adds a degree
of freedom. After eliminating the vector time component and the
lapse/shift constraints, the physical scalar velocity form is

    (s_dot+w*v_dot/Theta)^2/2
      + C*(sigma_dot+d*v_dot/Theta)^2
      + [J-(4/3)*q/h^2]*v_dot^2/Theta^2,
    C=4*zeta*q/(8+3*zeta*q)>0, h=(1+u^2)^3.

It has exactly one negative scalar direction whenever
q>3*J(u)*h(u)^2/4, at every finite u!=0 and positive normalized zeta.
Here J, Theta and w are the original coupled CD/M1 coefficients.
This is a direct configuration-velocity congruence in a regular
punctured chart, with the actual free chi retained. The threshold
tends to 3597/3200 at the crossing. A separate regular Hamiltonian
calculation there retains the time derivatives of the coefficients
and checks the baseline as a control; its instantaneous inertia is
not substituted for the physical punctured proof.

## What these rejections do not decide

An exact background is not sufficient for a healthy parent, and a
healthy-looking isolated Proca block is not its coupled spectrum.
Conversely, the positive-momentum condition is not an identified
low-frequency EFT band. The zero kinetic coefficient changes
constraint rank, and a small nonzero coefficient can move frequencies
outside a separately prescribed domain. No ghost-frequency cutoff,
loop remainder or controlled light-only EFT exclusion is claimed.

The physical curl coupling and the normalized one differ by
zeta=zeta_physical/(M^2*tau^2). A nonunit calibration checks the action
factor, vector normalization and isolated mass conversion. These
checks do not establish a heavy gap.

Neither the three literal additions nor published torsion-free
Minkowski spectra classify all unrestricted scalar-dependent
metric-affine parents. In particular, arbitrary kinetic contractions,
mixed connection traces and degenerate kinetic charts require their
own constraints. A bounded next task is to test constant mixtures
of the two independent projectively invariant distortion traces,
including their zero-Schur cases. It must not inherit a verdict by
continuing the nonsingular Proca formula through a rank change.

The adopted vacuum, finite-gravity and controlled matching V/G/B
obligations remain unresolved research. No intervention or new
assumption from the user is required to continue the local analysis.

## Verification

The frozen report covers 18 sources, 104 named exact identities
comprising 5529 scalar entries, 23 continuous/interface proof checks
and 24 rejected-input controls. All 142 ordinary tests pass in 86.48
seconds, and the separate read-only certificate CLI passes. Neither
uses the broad-regression GCD adapter; both random seeds are zero and
the diagnosed host faulthandler timer plugin is disabled.

The final adversarial action/domain review found two defects before
freezing: the full curvature kinetic matrix needed its trace projector,
and the formulation needed to distinguish the homothetic Maxwell
gauge from the other actions' projective quotient. Both were corrected
and replayed; neither changed the physical-sector signs. These are
symbolic certificates and written proofs, not Lean formalization or
external peer review.

The full P8 regression through this checkpoint passes **3945 tests in
686.61 seconds**. Only the not-yet-published S6.39 `kinetic/traces`
child is excluded from that checkpoint run. The existing opt-in exact
GCD runner passes all 128 untouched normalized-tuple comparisons and
reports 4996 domain fallbacks and 6509 exact descents. No polynomial,
assertion or frozen source was changed for this run.

The first broad run was interrupted after 1398.45 seconds with 1214
tests passed and no failures, while an older generic stationary-series
test remained in SymPy polynomial arithmetic. The successful retry
sets `PYTHONHASHSEED=0` at process start and resets SymPy's RNG to zero
before each test through an ephemeral pytest setup hook. It uses the
same exact runner and test assertions, with `--import-mode=importlib`
and `-p no:faulthandler`. This is a reproducibility/performance control,
not a mathematical change to the certificate.
