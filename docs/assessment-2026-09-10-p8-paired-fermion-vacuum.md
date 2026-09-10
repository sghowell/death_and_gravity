# P8: both paired fermion vacuum primitives

Original P8 remains OPEN. This follows the
[insertion local-reference checkpoint](assessment-2026-09-10-p8-finite-ms-insertion-references.md).
No user intervention is required for the current research.

## Both remaining primitive rows have assigned bounds

[S6.149](../problems/P8/s6/continuation/s6_149/FORMULATION.md)
evaluates the scalar and gauge fermion vacuum primitive rows,
including their complete assigned proper forests. The three
two-edge subcycles overlap; their counterterms are not multiplied
as though they were disjoint.

Dimension-symbolic vacuum tensors and both proper fermion
mass/kinetic counterterm insertions give finite massless-exchange
coefficients -19 NYm^4/Q^2 and +18 NaC_Fm^4/Q^2. Their full
scale dependence reproduces the independently calculated
quadratic mass anchors after fixed-scale differentiation.
Finite regulator/pole products are kept before MS subtraction.

The actual scalar exchange has mass one, not zero. Its entire
physical mass/residue subgraph is paired before the vacuum loop
is evaluated. Three exact beta anchors separate its finite
large-mass terms from a convergent remainder bounded by 12/T,
T=4m^2. The resulting scalar finite vacuum term is

    NY/Q^2 {-19m^4+T(13/4+pi^2/8)-47/18-pi^2/12+delta_V}.

The vector Ward identity makes the gauge-cycle counterterm
scaleless when its gauge legs are closed; the proper fermion
counterterms still contribute. Gauge exchange includes all
fourteen Dirac flavors, or 42 color states, not just the six
Yukawa-active states. The scalar contribution is bounded by
approximately 3.54193856196148e592 and the gauge contribution
by 7.41746690538194e593. Their sum is below 1e595 and below
1e-203 relative to the retained one-loop fermion vacuum reference.

All nine listed fermion primitive rows now have bounds in their
stated scopes. This is not a completed whole-model vacuum
counterterm or canonical amplitude. Assigned proper counterterms
are removed from the separate insertion ledger rather than
counted twice. Other scalar, source and matching terms remain.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
53 named identities, 53 scalar entries, 23 proof gates,
9 controls and 91 rejected inputs. Final private science passed
237 tests in 22.95 seconds; fresh repository science passed
237 in 23.86 seconds. Ordinary replay passed 262 tests in
2131.11 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 23180 tests in 3244.41 seconds.

All 551 captured test files were present and unchanged. Path-list SHA-256:

    3f3511a88ba36b4bd6eaa2765ce69edb4707206e44fe9fe91433cb714af891af

The full-run adapter passed 128 original tuple comparisons.
Counters were 34060 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.150. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 37036-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    32dcc39d645f6c9375280b99799516dce90c2982b41c694c210bf199793a3909

Independent tests include explicit Dirac matrices and Ward
identities, complex-regulator finite parts, fixed-scale mass
derivatives, the full dimensional scalar spectral identity,
Gamma duplication, three beta anchors, vacuum integral
representations and finite-remainder/pole-product checks.
All private cleanup preceded the final private run and freeze.
No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization or
independent peer review. Exact 22-file staging excludes later
continuations and unrelated P4/P9 changes.

## Remaining work

The isolated fixed-contact conversion and complete wineglass
interaction-forest MS conversion have passed their private,
repository and native checks; fresh independent regressions
are running. Double-bubble conversion is private research.
These are not counted as published results here.

Other scalar interaction families and full parameter/field
re-expansion, vacuum/source assembly including the first-source
square, complete canonical pole/error, physical higher-order
truncation, V contour/cut control, finite-gravity G and
common-parent B remain open. No nonperturbative gauge vacuum
or global effective-potential conclusion follows from these
fixed-order vacuum terms. No all-orders UV construction is
added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged. Original P8(b) and
original P8 remain OPEN.
