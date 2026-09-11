# P8: complete matched two-loop vacuum reference

Original P8 remains OPEN. This follows the
[complete matched Phi forward-coefficient checkpoint](assessment-2026-09-10-p8-complete-matched-gy14-two-loop-amplitude.md).
No user intervention is required for the current research.

## Complete vacuum ownership

[S6.158](../problems/P8/s6/continuation/s6_158/FORMULATION.md)
completes the two-loop vacuum reference in the same GY14
MS-interaction and physical-Phi reference prescription.

Before finite parts, the full regulated first stationary
source J1_D=-G T_D(1)/2 gives three reducible contributions:

    heavy-reducible dumbbell   -g T_D(1)^2/(8M),
    source-induced mass term  +g T_D(1)^2/(4M),
    first-source square       -g T_D(1)^2/(8M).

They sum to zero. Adding a source square on top of the
full-field 1PI vacuum would double count. Replacing the
regulated source by its finite part before squaring would
lose pole times positive-epsilon products.

The remaining complete scalar coefficient is the finite
Laurent part of

    -L T_D(1)^2/8 - g S_D(1,1,M)/4
    +g B_HL,D(-1)T_D(1)/2 + g T_D(M)/(4Q epsilon).

This includes the fixed physical-Phi mass/residue trace
and the proper heavy MS mass counterterm. The correlated
Phi slope insertion closes to a scaleless constant trace;
it is not dropped in isolation. Both complete fermionic
vacuum primitives of S6.149, including all their assigned
proper counterterms, are then added exactly once.
The gauge vacuum counts all 42 color/flavor states.

## Finite-part bounds and matching

Six ordered Schwinger sectors explicitly subtract the
sunset corner. On the regulator circle |epsilon|=1/16,
the compact remainder, gamma factors and massive powers
give

    |S_D(1,1,M)| < 4800 mF^(1/4)(3M)^(9/8)/Q^2.

The finite Laurent coefficient is bounded by the same
circle maximum. Exact rational power caps give a scalar
vacuum allowance approximately 8.62336e265, strictly below
1e275. Including both fermionic vacuum primitives gives

    complete second vacuum bound  7.77166076157809e593,
    ratio to first vacuum lower   3.15800818248253e-206.

Strictly the absolute second-order bound is below 1e595
and the stated ratio below 1e-203. Large dimensionful
vacuum coefficients are fixed renormalization references,
not a prediction of a nonzero physical vacuum energy.

The first Phi map changes interaction coefficients but
not the masses in the regulated first vacuum determinant.
Its constant Jacobian is scaleless. There is therefore
no omitted first-coordinate variation of that determinant,
and the second map acts on a zero tree vacuum reference.
The physical zero-energy condition fixes the finite
counterterm; no extra fitting freedom is introduced.

The complete second H-source is a separate next checkpoint.
Squaring its finite bound would not establish higher-order
regulated source-square references.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
50 named identities, 50 scalar entries, 25 proof gates,
9 controls and 147 rejected inputs. Final private science
passed 218 tests in 23.18 seconds; fresh repository science
passed 218 in 23.08 seconds. Ordinary replay passed 243
tests in 2168.04 seconds. Independent native CLI replay
passed. The full P8 snapshot passed 25067 tests in 3337.97
seconds with final exit code 0.

All 569 captured test files were present and unchanged.
Path-list SHA-256:

    69a370188981f8144a5411f30109f31d5962d809339f38845ddc4b91c7990529

The full-run adapter passed 128 original tuple comparisons.
Counters were 34087 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.159.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 41799-character native report was transferred losslessly
and all eighteen source hashes were verified. Report SHA-256:

    5867f46444e1227ea22cce2cfc70ac5fc4c8e28893e90ad5c2ae1e9632b22430

Independent tests check Gaussian source elimination,
regulated finite products on complex circles, all scalar
forest coefficients, six-sector analytic continuation and
rational mass-power bounds. Gaussian quadrature resolution
was increased during private development to meet the
unchanged accuracy tolerance. All cleanup preceded the
final private verification and immutable byte freeze.

No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The complete second stationary-H source has passed native
and direct science and is completing fresh regressions.
Full derivative-action transport with the matched physical
source is in private verification. Neither is counted
as published in this checkpoint.

Other required coordinate dictionaries, physical finite-EFT
truncation, V contours/cuts, finite-gravity G and common-parent
B remain open. A local vacuum reference is not a global
effective-potential or nonperturbative backreaction result.
No all-orders UV construction is added to the adopted
finite-EFT/necessary-positivity contract. Scoped P8(a) and
A.20-A.23 are unchanged; original P8(b) and P8 remain OPEN.
