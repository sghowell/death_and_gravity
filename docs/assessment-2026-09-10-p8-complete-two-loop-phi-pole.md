# P8: complete GY14 two-loop Phi normalization and local pole

Original P8 remains OPEN. This follows the
[regulated finite-field checkpoint](assessment-2026-09-10-p8-regulated-finite-field-covariance.md).
No user intervention is required for the current research.

## Complete quadratic ownership and regulator-complete field map

[S6.156](../problems/P8/s6/continuation/s6_156/FORMULATION.md)
assembles all thirty-two scalar refinements and all three
fermionic quadratic families with their assigned proper
counterterms. Physical inner OS conditions and direct-MS
fermion references are retained.

The full first coordinate shift is q1,D=-k_D Ephi q,
not just its finite coefficient. In the bare field
identity the positive-epsilon shift of the first MS pole
cancels the matching field-pole product. The resulting
second normalization is

    t_MS=t_H-k0 r0.

There is no residual z11 k1 term. Both canceled terms
remain displayed in the calculation and negative controls.

Re-expressing the first interaction poles at the full
regulated shifted coordinates gives

    G2=(k0^2-t_MS)G+k1 L G/Q,
    L2=(3k0^2-2t_MS)L+3k1 L^2/Q,
    M2=k1 g/Q.

This includes the induced heavy-mass shift and the
fundamental G1 square when forming g2. It completes
the bare scalar map that S6.155's fixed-input coefficient
identity did not by itself supply.

The complete second normalization is bounded by
approximately 3.34630300149847e-19, below 1e-18.
Its corresponding second tree-b2 allowance is approximately
6.69260600299694e-19 and strictly below 1e-18.

## Complete canonical unit-disc light pole

The direct fermion primitive estimate is extended from
its integrated radius-two holomorphic bound to the
unit disc by an explicit Cauchy-tail estimate. The scalar
and inserted-fermion OS families, and both first-parameter
variations, are retained in the same assembly.

For |s-1|<=1 the complete canonical inverse through two
loops has

    D2(s)=(s-1)+Pi_R(s),
    |Pi_R(s)|<1e-18 |s-1|^2.

Its only pole on this disc is therefore the mass-one,
unit-residue light pole. The formal field normalization
also stays positive. This is a fixed-order local result,
not a global or exact all-orders spectral theorem.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
51 named identities, 51 scalar entries, 24 proof gates,
9 controls and 96 rejected inputs. Final private science
passed 173 tests in 21.96 seconds; fresh repository science
passed 173 in 21.96 seconds. Ordinary replay passed 198
tests in 2199.85 seconds. Independent native CLI replay
passed. The complete P8 snapshot passed 24644 tests in
3298.13 seconds with final exit code 0.

All 565 captured test files were present and unchanged.
Path-list SHA-256:

    24c79abba5f66e75d6456a7ef275f21280dd74e9145675ab70c4d07650c9599b

The full-run adapter passed 128 original tuple comparisons.
Counters were 34084 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.157.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 65429-character native report was transferred losslessly,
and its eighteen source hashes were independently verified.
Report SHA-256:

    a75a1b9dc9f16e98e9d58500b3b36cb48f17d2744f3973c95108234facf37080

Independent tests use literal two-regulator field ratios,
full bare-vertex substitutions with first-pole re-expression,
complex-circle forward extraction including G1^2 and M2,
and analytic-tail/pole controls. The full regulated
coordinate correction was made and checked in the private
draft before freeze. No frozen bytes were changed.

Written analytic arguments and exact replay are not
formalization or independent peer review. Exact 22-file
staging excludes later continuations and unrelated P4/P9.

## Remaining work

The matched two-loop forward-coefficient and complete
vacuum-reference/source-square checkpoints have passed
direct science and native checks and are completing
fresh regressions. Neither is counted as published here.
The remaining complete second H-source assembly is in
private development.

Other required coordinate dictionaries, physical finite-EFT
truncation, V contours/cuts, finite-gravity G and common-parent
B remain open. No all-orders UV construction is added to the
adopted finite-EFT/necessary-positivity contract. Scoped P8(a)
and A.20-A.23 are unchanged; original P8(b) and P8 remain OPEN.
