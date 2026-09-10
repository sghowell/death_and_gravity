# P8: both scalar and gauge quartic chord primitives bounded

Original P8 remains OPEN. This follows the
[opposite-chord bound](assessment-2026-09-10-p8-fermion-opposite-chord.md).
No user intervention is needed for the current research.

## Complete paired vertex subsets

[S6.142](../problems/P8/s6/continuation/s6_142/FORMULATION.md)
bounds the remaining twenty-four vertex-chord words in each quartic
fermion primitive. Together with the disjoint self-energy and opposite
subsets, all sixty words per primitive are now covered.

Each vertex word has one ultraviolet-divergent proper Yukawa vertex.
Its full MS kernel is the finite local anchor plus the convergent
difference from its zero-external-momentum integrand. Scalar and gauge
signs and counterterms are retained. In particular, finite dimensional
numerator contributions in the local anchor are not dropped.

The remaining integration is divided into low and high inner momentum
relative to the outer massive norm. The low region controls both raw
and subtracted kernels. In the high region only their paired difference
is integrated. Resolvent identities give the needed extra decay.
No separately divergent high-momentum piece is assigned a finite bound.

The fourth-and-higher soft Taylor tail is projected before the
integrations. Complete permutation symmetry makes the lower-degree
terms irrelevant to b2. Exact radial and logarithmic integrals give

    E_vertex <=5e14 N Y_hi^2(Y_hi+4 a_hi C_F)/(Q_lo^2 mF^4).

At the shared reference, E_vertex is approximately
5.27382729154e-1404. Adding all three disjoint chord subsets gives
1.26571854997016e-1403, or 3.16429637492539e-804 relative to the tree.
This bounds the complete scalar_Phi4_W0_F4 and gauge_Phi4 primitive
rows, not the complete matched two-loop four-point amplitude.

## Scope and independent verification

Five other primitive rows are wholly unevaluated at this checkpoint.
The two older paired sectors are unchanged. Their reference conversions,
other parameter/field contributions, and the complete error budget
are not inferred from these small primitive bounds.

The report pins 18 source/proof/test files and 20 fields: 66 named
identities, 81 scalar entries, 29 proof gates, 9 controls and
66 rejected inputs. Final private science passed 282 tests in
21.53 seconds; fresh repository science passed 282 in 22.82 seconds.
Ordinary replay passed 307 tests in 2613.68 seconds.
Independent native command-line replay passed.
The complete P8 snapshot passed 21458 tests in 3423.48 seconds.

All 537 captured test files were present and unchanged. Path-list SHA-256:

    2f45d98a6835421b861c6fbdb24d4546a3c208c5a482d599c9e0d2f1b6a366aa

The adapter passed 128 original tuple comparisons. Full-run counters:
34033 domain fallbacks, 7340 exact descents and 88 mixed fallbacks.
The full snapshot predates S6.143. Only full regression uses the
audited exact GCD adapter. Native, ordinary, command-line and direct
science retain unmodified SymPy, with interpreter-only allowances.

The 57580-character native report was transferred losslessly and all
eighteen source hashes checked independently. Report SHA-256:

    2b0f5a84ac0e2fc1a0a5ffaba7c4930c4a0006c9d5d7cc66853475b21db33ea0

Tests independently enumerate proper graph cycles, all labelled words,
complex matrix resolvent bounds, low/high radial moments and frontier
mutations. Analytic proofs are written arguments, not formalization
or independent peer review. Exact 22-file staging excludes later
continuations and unrelated P4/P9 work.

## Continuing work and nonclosure

The separate finite outer MS-reference conversion has passed direct
science and native replay; fresh independent replays are running.
A private derivation of the remaining mixed quartic primitive now has
passing initial algebraic residuals and conservative reference bounds,
but has not been frozen or independently replayed.

The remaining primitive vacuum/quadratic rows, full matching, canonical
pole and field/parameter conversions, total two-loop and truncation
errors remain open. So do V contours, finite-gravity G and common-parent B.
No all-orders UV construction is added to the adopted finite-EFT and
necessary-positivity task. Scoped P8(a) and A.20-A.23 remain unchanged.
Original P8(b) and original P8 remain OPEN.
