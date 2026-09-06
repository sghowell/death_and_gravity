# P8 continuation audit after A.12 and S6.7.COMPOSITE

Recorded 2026-09-06. This extends the
[A.11/S6.6 assessment](assessment-2026-09-06-p8-prepared-see-composite-bounce.md)
without changing the original P8 statement or immutable certificate scopes.
Original P8 remains open.

## Two substantive screens

| Track | Established | Remaining boundary |
|---|---|---|
| [A.12](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei/FORMULATION.md) | A fresh absolute all-Hadamard QSEI on the actual source-free A.11 solution, using only quantitative C1 potential bounds | The certified short interval provably cannot meet the comoving focusing criterion; longer-domain and realistic-field results remain open |
| [S6.7.COMPOSITE](../problems/P8/s6/matching/composite/perturbations/FORMULATION.md) | Literal source-aware tensor/vector actions, finite vector constraints at zero tensor stiffness, and negative vector principal signs for both frozen bounces | An independent parent scale reverses the leading sign in a different candidate; scalar, cutoff, matching and positivity checks remain open |

Certificate SHA256 values:

    A.12: f72698170ee7ca0468b42ceaa098e69e1a26054c6ffdbf6321f1968eaefd6584
    S6.7.COMPOSITE: 018295c17af163cc80801e3dc2a9d7309184e4d1f5eca319b2a59ea85f064090

A.12 pins A.11, hashes 18 files and checks 23 symbolic identities plus an
independent Fraction reconstruction. S6.7 pins S6.6, hashes 17 files and
checks 49 symbolic and eight independent Fraction coefficient identities.
The proofs include analytic arguments not covered by symbolic residuals;
neither gate is Lean-formalized. All previous hashed subtrees are unchanged.

Validation: **1,301 P8 tests passed** in 405.52 seconds, including all 54
A.12 and 37 S6.7 tests. Both read-only certificate CLI replays and the full
P8 Ruff check pass. P4/P9 work is excluded.

## What the exact-solution QSEI does and does not buy

A.12 does not assume that A.9's higher-jet estimate survives the metric
change. A new decomposition separates a real Fourier history of the first
potential derivative from the complex mode remainder. Momentum Parseval
controls that history; full sampling-frequency Parseval controls complex
products. The distinction matters: halving a complex transform's spectrum
would be invalid. No unproved uniform second or third potential derivative
cap is imported from the old metric.

The actual reference state's effective energy is positive by the exact
Einstein equation minus the fixed ordinary radiation and the A.11 metric
distance. This turns the difference bound into an absolute bound for every
Hadamard target state, not only homogeneous or nearby states. The sufficient
coefficient is `2*hbar/(16*pi^2)`; its derived rational precursor is closer
to the flat coefficient, but no optimum is claimed.

The free half-slab has proper span at most `(3L/2)*T0`, `L=10^-10`.
An independent actual-curvature index estimate exceeds the available
comoving contraction for every endpoint test function wholly within this
span. Thus failure to apply the focusing theorem is not cured merely by
optimizing the QSEI coefficient there. It also does not imply completeness
or exclude other hypersurfaces, congruences or longer domains.

The next useful calculation is the full linear response, retaining the
stiff Einstein term as well as the logarithmic quantum response, to see
whether a resummed continuation estimate is viable. Iterating a short-slab
existence theorem without accumulated bounds is not such an estimate.

## What the vector screen changes

The relative tensor stiffness includes the composite matter pressure.
The shift stiffness includes variation of both the effective volume and
inverse lapse. Exact constraint elimination uses the finite algebraic
shift coefficient Xi; the shorthand `mu/cV^2` is not divided at its zero.
The physical clock and ruler give the principal speed `(r/s)^2*mu/Xi`.
The literal unreduced-action calculation also records a bounded
normalization discrepancy in the source's printed reduced vector formula.

For the two frozen S6.6 candidates, exact negative quadratic bounce jets
and analyticity prove a punctured interval of negative vector principal
coefficient, with positive vector inertia. This is sufficient to reject
their formal strict vector-principal-health requirement. It is not a
calculation of large amplification below an established physical cutoff.
Nor does vanishing algebraic tensor stiffness alone determine a canonical
rolling mass gap; normalization derivatives have been kept explicitly.

The earlier identification of parent mass `m` with inverse duration
`1/tau` belongs to those frozen examples. Independently varying `m*tau`
reverses the leading CD vector sign above `sqrt(24)`. The value 5 supplies
an exact positive control, not a healthy-bounce theorem. A new candidate
must independently satisfy its background domain, regular scalar
constraints, gradients and matching conditions. Published high-frequency
kinetic formulas are not a substitute for that full analysis.

## Closure decision and next work

The original realistic-field cosmological focusing and controlled
matter/positivity obligations remain. No missing user choice, credential,
external approval or compute budget has been identified. The next work is
the long-domain response calculation on A and the independently scaled
background plus regular scalar analysis on B. These are mathematical
research tasks, not grounds for marking P8 finished or closed.
