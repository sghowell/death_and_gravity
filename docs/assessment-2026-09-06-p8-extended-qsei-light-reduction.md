# P8 after actual SEE extension, fresh QSEI and rolling tensor reduction

Recorded 2026-09-06, following the
[resummed-response/cone checkpoint](assessment-2026-09-06-p8-resummed-response-cone-boundary.md).
Original P8 remains open. These results close three specific proof gates,
not the realistic-field cosmological theorem or controlled positivity
matching required by the original problem.

## Verified new scope

| Gate | New result | Remaining boundary |
|---|---|---|
| [A.14](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/extension/FORMULATION.md) | Actual smooth full SEE solution extended from L0=10^-10 to L=10^-6, preserving the identical source, state, data and original source-off time | Local, not macroscopic continuation or a new physical no-runaway prescription |
| [A.15](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/extension/qsei/FORMULATION.md) | Fresh absolute coefficient2 QSEI for all Hadamard targets and proper H2 samplers on the longer source-free domain, with a signed actual-reference bound | The actual comoving index trigger still fails; realistic-field cosmological incompleteness is not established |
| [S6.9.COMPOSITE](../problems/P8/s6/matching/composite/perturbations/cones/matching/FORMULATION.md) | Exact time-dependent two-tensor action/retarded reduction and all-duration hierarchy screens at two specified initial states | No general-state, finite-band nonadiabatic, scalar, cutoff or UV matching verdict |

Certificate SHA256 values:

    A.14: 93c0ac119240a8c835da52f456d6e7cdc846e872728717de869e8fb823275ec1
    A.15: 171e600f90a9c5e0240addaceaa647260001fe4f11bc4c9045d4945193de069e
    S6.9.COMPOSITE: 36c899d2b82957f2fc83df8e7bc66582cbf77829fbed8b64724d2d39c2086e48

A.14 hashes 18 files and checks 11 exact identities; A.15 hashes 14 files
and checks 25 pinned-general/new identities; S6.9 hashes 14 files and
checks 33 identities plus four separate Fraction coefficient identities.
Every new numerical gate is recomputed with exact rationals, with separate
physical/covariant audits (9, 8 and 9 tests respectively). Functional
analysis, state existence and continuous-time theorem implications remain
written proofs, not proof-assistant formalizations or finite sampling.

Validation: **1,498 P8 tests passed in 393.93 seconds**, using
`PYTHONHASHSEED=0` and diagnostic `faulthandler_timeout=60`. This includes
45 A.14, 37 A.15 and 39 S6.9 tests. All read-only certificate replays,
full P8 Ruff, local Markdown artifact links and diff checks pass. Two
diagnostic slow-test stack dumps came from unchanged S5 symbolic algebra;
neither was a test failure. Immutable ancestors and unrelated P4/P9 work
are unchanged. Unfinished follow-on children are outside this checkpoint.

## A: a larger interval, with the same physical solution

The complete resummed map uses the causal inverse with exponential weight
sigma=2*10^6. Its weighted contraction is below 3*10^-6; the actual
fixed-point distance is below 9*10^-8. The exp(2)<9 conversion is included
in every required pointwise barrier, and the auxiliary bounds are proved
before being used. Smoothness retains the same flat-start neighborhood
and original history. The actual Hadamard state, full stress and density
constraint follow from the smooth fixed point, and old uniqueness proves
agreement on the original slab.

The source still turns off at the original L0/2=5*10^-11. Thus the new
source-free conformal interval is almost the whole L=10^-6 extension,
not just its new final half. No reset of quantum covariance or energy,
rescaling of the source cutoff, or deletion of a growing causal pole is
used. The weight does not justify arbitrary longer intervals.

The new QSEI uses the actual extended geometry and state. Its C1-only
mode lemma retains the full preparation history, while the sampler
norms are recalibrated for the new length. The old positive reference
credit cannot be transferred. Instead, the weighted metric difference
and the actual unforced SEE give
`E_reference>=-hbar/(40pi^2 T0^4)`. Proper Poincare bounds add only
`(2/5)L^4` to the derivative coefficient, leaving the full absolute
coefficient strictly below the rounded value2. Target states need not
be homogeneous, quasifree or solutions of the SEE.

The sufficient Q2/available-duration^2 ratio improves from 160000000
to at least2/5. Nevertheless, the actual-curvature index form on this
interval is greater than the available initial contraction, in either
time orientation. That result cannot be repaired by merely rounding
the QSEI coefficient more tightly. There is no cosmological focusing,
maximal-extension or incompleteness conclusion from this local interval.

## B: algebraic mass is not the rolling heavy response

The exact physical-clock action transforms to weighted and relative
tensor fields. The canonical transformation retains its time boundary,
moving eigenvector, normalization curvatures and gyroscopic coupling.
The heavy equation gives `H=H_hom-G_R B l`; substitution in the light
equation retains `B*H_hom`. A retarded kernel is not varied as if it
were an ordinary symmetric single-copy action. The first inverse-mass
representative has an explicit differential residual, whose conversion
to a physical error would require a propagator and initial-data bound.

At the beta2-model CD center with rho0=1/2, y0=1 has zero algebraic
relative mass. At y0=2 the algebraic mass is positive and the leading
locked cone is subluminal, but the ratio of either physical curvature
or eigenvector mixing to that mass squared exceeds50/9 for every
nonnegative CD acceleration A. Changing m*tau alone therefore cannot
produce the specified slow mass-led hierarchy at these data.

An exact fixed-charge zero-momentum Routh reduction supplies a needed
countercontrol: its frequency can be positive even at zero algebraic
mass. A negative diagonal canonical potential or zero stiffness alone
is not a physical growth or finite-band spectral-gap theorem. None of
the specified-data tests excludes every other state or a controlled
nonadiabatic reduction, and none integrates the scalar/vector sectors.

## Next work and closure boundary

The original problem asks for realistic-field cosmological QEIs and a
stable-bounce classification with controlled matter/UV-positivity
applicability. A nearby exact minimal-scalar SEE solution is a useful
chosen route, not a logical prerequisite for every possible singularity
theorem. A separate primary-source photon-field audit is therefore being
considered alongside further actual-solution control; it does not yet
establish a Maxwell or interacting-field theorem.

On the matching track, large-asymmetry states are a concrete remaining
test rather than an assumed extension of the y=1,2 screens. Their full
normalization, state/source projection and finite-time dynamics need to
be checked before drawing a hierarchy or instability conclusion.
Any such follow-on work is separate from this frozen checkpoint.

No missing credential, external approval or unresolved user choice
currently blocks those in-scope calculations. Scientific obligations
remain; P8 is not finished or closed.
