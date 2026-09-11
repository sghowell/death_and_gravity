# P8: complete matched GY14 two-loop Phi forward coefficient

Original P8 remains OPEN. This follows the
[complete Phi normalization/pole checkpoint](assessment-2026-09-10-p8-complete-two-loop-phi-pole.md).
No user intervention is required for the current research.

## Complete four-point ownership and matching

[S6.157](../problems/P8/s6/continuation/s6_157/FORMULATION.md)
assembles all 192 scalar MS-interaction refinements and
all four fermionic quartic families. Each proper and
overall reference is owned once; the full heavy kernels,
physical inner OS conditions and direct-MS fermion
counterterms remain in their assigned families.

Use the hybrid representation with physical Phi OS
conditions and MS interaction poles. The finite scalar
one-loop functional has Phi-coordinate weight four;
the fermionic box has weight two. The first coordinate
variation is consequently

    -k0[4 A_H1,scalar+2 A_H1,fermion].

The complete renormalized first-loop sum is finite
before this variation. This does not permit dropping
positive-epsilon terms in the separate bare map:
S6.156's full G/L/M conversion, including k1, M2 and
the fundamental G1 square, remains in the second tree
contribution. No extra raw-MS field or LSZ factor is
applied to an already hybrid/canonical amplitude.

The isolated sigma cancellation is a pointwise
common-regulator functional identity. Composing it
with the field map preserves its order without assuming
that the two parameter directions commute. A separate
sigma insertion is therefore not added again.

## Complete through-two-loop coefficient

For T=2g/(M-2)^3=4lambda, the complete second-order
coefficient is

    b2^(2)=b2_H,raw^(2)
      -k0[4 b2_H,scalar^(1)+2 b2_H,fermion^(1)]
      +T[3k0^2-2t_MS+(k1/Q)(2L-3g/(M-2))].

The four raw scalar families and four fermion rows,
the finite first-loop variation and the second tree
map exhaust the named contribution. Local Phi2/source
references enter through the already fixed physical
pole; Phi parity forbids an omitted internal light
tree exchange.

The exact rational allowances, printed approximately,
are

    raw second relative       6.73899810711317e-8,
    first-coordinate cross    6.87384305853845e-213,
    second tree-field map     6.69260600299694e-19,
    complete second relative  6.73899810718010e-8,
    first plus second         5.14424816887231e-7.

Strictly, the complete second-order relative bound is
below 1e-7 and the combined bound below 1e-6. The formal
two-loop b2 polynomial is positive for 0<=h<=1. The
complete canonical unit-disc light pole of S6.156 is
unchanged.

This bounds the computed perturbative coefficients,
not the physical omitted higher-loop remainder.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
41 named identities, 41 scalar entries, 25 proof gates,
9 controls and 90 rejected inputs. Final private science
passed 155 tests in 21.17 seconds; fresh repository science
passed 155 in 20.91 seconds. Ordinary replay passed 180
tests in 2239.79 seconds. Independent native CLI replay
passed. The complete P8 snapshot passed 24824 tests in
3362.16 seconds with final exit code 0.

All 567 captured test files were present and unchanged.
Path-list SHA-256:

    62c502d3d607797119e93ff28e5d526410cd783b525b1ec0066506700582d0a1

The full-run adapter passed 128 original tuple comparisons.
Counters were 34102 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.158.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 106382-character native report was transferred
losslessly and all eighteen source hashes were verified.
Report SHA-256:

    ca620649aa67631a60681b9e32f2cac204cae0f530230261a9490d1064b5651f

Independent tests extract the forward coefficient from
the literal fully mapped amplitude on a complex contour,
retain the induced heavy-mass shift and cubic square,
check noncommuting coordinate composition, and detect
an extra LSZ factor or an early regulator projection.
All cleanup preceded final private verification and freeze.

No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The complete vacuum reference and complete second H-source
checkpoints have passed direct science and native checks
and are completing fresh regressions. Neither is counted
as published here. Transport through the full derivative
action with its matched physical source is under review.

Ordinary off-shell coordinate normalization is not assumed
uniform at arbitrarily high momentum. Other required
coordinate dictionaries, physical finite-EFT truncation,
V contours/cuts, finite-gravity G and common-parent B
remain open. No all-orders UV construction is added to the
adopted finite-EFT/necessary-positivity contract. Scoped
P8(a) and A.20-A.23 are unchanged; original P8(b) and P8
remain OPEN.
