# P8: complete wineglass interaction-forest MS conversion

Original P8 remains OPEN. This follows the
[fixed-contact conversion checkpoint](assessment-2026-09-10-p8-fixed-contact-ms-conversion.md).
No user intervention is required for the current research.

## All sixteen wineglass refinements

[S6.151](../problems/P8/s6/continuation/s6_151/FORMULATION.md)
converts the entire sixteen-refinement wineglass interaction
forest to MS at mu=mF, at fixed canonical reference coordinates.
The local anchor is derived from the full dimensional equal-mass
sunset derivative. Its sector-subtracted parameter integral fixes
a finite constant j, bounded by |j|<=22.

The proper-MS zero-external reference has finite part

    [j+ell+ell^2/2+pi^2/12]/Q^2,

where Q=16 pi^2 and ell=log(mu^2). This differs from the old
recursive reference. The epsilon coefficient of the one-loop
reference must be retained in its pole product.

An anchor alone does not complete this conversion. The finite
inner-scale change is integrated over the full outer light/heavy
propagators, including both heavy numerator terms and every
channel. Its b2 bound is approximately 1.66616720420620e-611.
Adding the old family and exact local heavy-tree variation gives
a complete wineglass-MS bound of 2.50439946558600e-611,
or 6.26099866396501e-12 relative to tree.

The equivalent linear second-order G and L reference shifts
are fixed by that same local amplitude term. They are an
alternative representation, not an additional contribution.
No extra LSZ factor or fitted finite contact is introduced.

This is a complete assigned interaction-family conversion,
not the global field/coupling map or its re-expansion terms.
The nine fermionic primitive rows remain unchanged.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
48 named identities, 48 scalar entries, 23 proof gates,
9 controls and 88 rejected inputs. Final private science passed
183 tests in 28.52 seconds; fresh repository science passed
183 in 28.26 seconds. Ordinary replay passed 208 tests in
2053.72 seconds. Independent native command-line replay passed.
The complete P8 snapshot passed 23572 tests in 3378.63 seconds.

All 555 captured test files were present and unchanged. Path-list SHA-256:

    fa7998dc294e388ec38d418be58d4fa62e7219814faed831efe0c7bb92eae5a9

The full-run adapter passed 128 original tuple comparisons.
Counters were 34074 domain fallbacks, 7340 exact descents and
88 mixed fallbacks. This snapshot predates S6.152. Only full
regression uses the separately audited exact GCD adapter;
native, ordinary, CLI and direct science use unmodified SymPy
with interpreter-only allowances.

The 31320-character native report was transferred losslessly,
and all eighteen source hashes were independently verified.
Report SHA-256:

    347356124e5cda2c5d2c68510c5c2f7325414f0559b49219c068b5d8e81295ad

Independent tests compare the sector representation with a
separate bubble/radial hypergeometric integral at nonzero
regulator. Complex-circle extraction checks the regulated
double-pole, single-pole and finite terms at three scales.
Other tests retain both heavy terms, exact tree-coordinate
variation and the full three-channel forward coefficient.

All private cleanup preceded final private verification and
freeze. No frozen scientific source or report changed.
Written analytic proofs and exact tests are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The double-bubble and scalar insertion conversions have
passed direct/native checks and are completing their fresh
regressions. The scalar quadratic finite-MS-slope continuation
is in private development. None is counted as published here.

Global field/coupling matching and cross terms, full GY14
vacuum/source and canonical amplitude/pole assembly, physical
higher-order truncation, V contour/cut control, finite-gravity
G and common-parent B remain open. No all-orders UV construction
is added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged. Original P8(b) and
original P8 remain OPEN.
