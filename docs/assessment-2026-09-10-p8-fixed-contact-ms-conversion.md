# P8: isolated fixed-contact MS conversion

Original P8 remains OPEN. This follows the
[paired fermion vacuum checkpoint](assessment-2026-09-10-p8-paired-fermion-vacuum.md).
No user intervention is required for the current research.

## One exact nonlocal matching cancellation

[S6.150](../problems/P8/s6/continuation/s6_150/FORMULATION.md)
discharges only the -sigma direction of the previously fixed
finite quartic contact's MS conversion. The scale, other
couplings and canonical field convention are held fixed.

The actual contact is affine, sigma_D=-k_D L+c_D, with a
common holomorphic dimensional lift. Its exact isolated
coordinate inverse is

    Lstar+h sigma_D(Lstar)=L,
    Lstar=(L-h c_D)/(1-h k_D).

The order-two inverse term is retained. The inherited
second-order local references are sigma_D times the L
derivative of the first-order references. They cancel the
corresponding first-order re-expansion before finite parts
are taken. Truncating sigma_D in only one pole-containing
leg would leave a finite defect.

The same identity cancels the complete assigned nonlocal
contact insertion against the re-expanded one-loop amplitude,
including both heavy triangles and all three channels.
This is not inferred from the zero b2 of a bare contact:
the contact's loop insertion is independently nonzero.

In this isolated direction the first field correction and
heavy source are unchanged. The quadratic insertion is a
local tadpole fixed by the existing physical mass condition.
No additional noncontact G2 or M2 shift, tunable finite
reference or external LSZ factor is introduced.

The inverse is nonsingular at the actual parameters. The
first and second coordinate-shift bounds are approximately
8.73114913702011e-409 and 1.66533453693773e-612; the remaining
coordinate tail is below 4e-816. That last estimate is an
algebraic inverse tail, not a physical higher-loop bound.
Only one isolated matching obligation advances; all nine
primitive rows remain unchanged.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
66 named identities, 66 scalar entries, 23 proof gates,
9 controls and 49 rejected inputs. Final private science passed
159 tests in 11.82 seconds; fresh repository science passed
159 in 11.89 seconds. Ordinary replay passed 184 tests in
2188.79 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 23364 tests in 3341.30 seconds.

All 553 captured test files were present and unchanged. Path-list SHA-256:

    72a37aa572bddbfde40fdf94c5d14cad2307d68ab7c2bd7f2e0c82f239a62459

The full-run adapter passed 128 original tuple comparisons.
Counters were 34083 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.151. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 35999-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    4f3154ca66df6c4fd39426bd87081f1f15078e687822908c2e718f80ca00e673

Independent tests include radial and Feynman-parameter
integrals at nonzero regulator, nonconstant affine maps,
literal three-channel amplitude/coupling expansions, both
heavy-term mutations and epsilon-times-pole finite defects.
An initial private test exposed a Fraction-input helper
mismatch; the checkpoint uses the project's Fraction-capable
exact validator. All fixes preceded the final private run
and freeze. No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The wineglass and double-bubble interaction-forest conversions
have passed private, repository and native checks; fresh
regressions are running. The scalar OS insertion conversion
has passed its final private checks and is being frozen.
None is counted as published here.

The remaining global field/coupling map and cross terms,
full GY14 vacuum/source and canonical amplitude/pole assembly,
physical higher-order truncation, V contours/cuts, finite-gravity
G and common-parent B remain open. No all-orders UV construction
is added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged. Original P8(b) and
original P8 remain OPEN.
