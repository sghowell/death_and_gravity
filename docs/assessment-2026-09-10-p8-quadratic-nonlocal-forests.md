# P8: both quadratic primitive nonlocal remainders are bounded

Original P8 remains OPEN. This follows the
[complete mixed quartic checkpoint](assessment-2026-09-10-p8-complete-mixed-fermion-quartic.md).
No user intervention is needed for the current research.

## What advances

[S6.145](../problems/P8/s6/continuation/s6_145/FORMULATION.md)
bounds the complete nonlocal on-shell remainders of
scalar_Phi2_W0_F2 and gauge_Phi2. Each sector contains three cyclic
words: two self-energy words and one vertex word. The full
proper-forest ledger includes both overlapping vertex subgraphs
and the whole fermion cycle. Its local mass term is removed only
by the explicitly stated on-shell affine projection.

The complete previously derived MS self-energy kernels are used.
For a vertex word, both overlapping proper subtractions are paired
with the graph, and both finite MS anchors are restored. There is
no product counterterm for overlapping subgraphs.

The joint momentum integral is divided into a comparable-momentum
region and both unequal-momentum regions. Each complementary
subtraction retains its own valid soft radius; neither a hard
subgraph's pointwise estimate nor a single-momentum bound is
silently extended over the whole integral. The sum is integrable
at the massless chord diagonal and at both ultraviolet ends.

For |s-1|<=1/2, the absolute on-shell divided remainder is bounded
by approximately 3.27434727200701e-800, below 1e-799. The bound
follows from the complete soft tail on |s-1|<=2, including its
affine on-shell projection. It is not a bound on the subtracted
finite mass or slope references.

Only the two nonlocal quadratic rows advance. The four quartic
rows and the older paired quadratic sector remain unchanged.
Two vacuum rows are still wholly unevaluated.

## Verification

The report pins 18 source/proof/test files and 20 fields:
68 named identities, 74 scalar entries, 24 proof gates, 9 controls
and 64 rejected inputs. Final private science passed 215 tests in
21.96 seconds; fresh repository science passed 215 in 21.98 seconds.
Ordinary replay passed 240 tests in 2525.74 seconds.
Independent native command-line replay passed.
The complete P8 snapshot passed 22210 tests in 3194.34 seconds.

All 543 captured test files were present and unchanged. Path-list SHA-256:

    e3f06cb382a5d670dda2774fc67263ec7d6cdc6f2e08e7cd2f17eff0383de5dc

The full-run adapter passed 128 original tuple comparisons, with
34019 domain fallbacks, 7340 exact descents and 88 mixed fallbacks.
This snapshot predates S6.146. Only full regression uses the
audited exact GCD adapter. Native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 30464-character native report was transferred losslessly.
All eighteen source hashes were independently checked. Report SHA-256:

    508c59cf4c3a58222bd3afd5af28b242643b6a5fe004e869af335b59e8a4d20f

Independent checks include cyclic enumeration, proper-cycle
overlap, exact radial antiderivatives, numerical integrals,
Cauchy projection and frontier mutations. One incorrect private
test equating distinct radial pieces was replaced by their
actual integral formulas before the final private run and freeze.
No frozen scientific source or report was changed.
The proofs are written arguments with exact algebra, not
formalization or independent peer review. Exact 22-file staging
excludes later continuations and unrelated P4/P9 changes.

## Remaining work

The next slope checkpoint has passed private, repository and
native checks; its fresh independent replays are running. The
finite mass calculation is still private and is not counted here
as an independently certified result.

Finite local references, vacuum rows, the other counterterm and
parameter conversions, complete canonical matching and the full
two-loop pole/error remain open. Higher-loop truncation control,
V contours, finite-gravity G and common-parent B remain open.
No all-orders UV construction is added to the adopted
finite-EFT/necessary-positivity contract. Scoped P8(a) and
A.20-A.23 are unchanged. Original P8(b) and original P8 remain OPEN.
