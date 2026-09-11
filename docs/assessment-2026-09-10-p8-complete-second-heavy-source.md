# P8: complete second H-source and fixed-order reference assembly

Original P8 remains OPEN. This follows the
[complete two-loop vacuum reference](assessment-2026-09-10-p8-complete-two-loop-vacuum-reference.md).
No user intervention is required for the current research.

## Complete stationary-source ownership

[S6.159](../problems/P8/s6/continuation/s6_159/FORMULATION.md)
completes the second stationary-H source in the unchanged
GY14 MS-interaction and physical-Phi reference prescription.

At fixed physical Phi mass one, the scalar covariance
insertion integrates to

    Tinsert_scalar,D=-g S_a/2-g B_HL,D(-1) I2+r_D T_D(1),
    S_a=partial_a S_D(a,a,M)|_(a=1),
    I2=-partial_a T_D(a)|_(a=1).

Both light lines contribute to S_a. Independently,
differentiating the vacuum functional with respect to H
gives the same source when its physical reference
counterterms are held fixed. Retuning the mass or residue
conditions at each H would change the theory and its
one-point function. The local L terms cancel in this
fixed-reference derivative.

The complete second source includes this scalar insertion,
the full fermion covariance insertion of S6.148, and the
proper cubic MS counterterm tadpole:

    J2,H=Fin[-G(Tinsert_scalar,D+Tinsert_fermion,D)/2
            -deltaG1,D T_D(1)/2].

The proper cubic term is not already inside the covariance
insertion. Gauge vacuum pieces and the remaining fixed
fermion/heavy reference terms are H-independent at Phi=0.
All source/reducible cancellations use the full regulated
first source, as established in S6.158.

The full first Phi-coordinate re-expression of the first
source adds

    -G[k0(ell+1)+k1]/(2Q).

The k1 term multiplies its first pole and is retained.
No second field factor acts on a nonzero tree source:
J0=0, and the bare H/source coordinate is unchanged.

## Bounds and assembled references

Six differentiated Schwinger sectors and the split
r0=(3M)^(-1) give the explicit regulator-circle bound

    |S_a| <28000 mF^(1/4)(3M)^(3/16)/Q^2.

The split controls the endpoint without introducing an
unnecessary full power of M. Rational power caps bound
all remaining scalar/source groups. Scalar plus cubic
terms are below 1e80; the first-coordinate source term
is below 1e-209. The complete source is dominated by
the already assigned fermionic insertion:

    complete second H-source upper  4.77877709968590442e188,
    upper divided by M              2.44673387503918306e-9,
    induced Phi mass reference      2.98673568730369026e-13.

Strictly these are below 1e189, 1e-8 and 1e-12 respectively.
They are fixed counterterm/reference allowances; physical
H remains at its stipulated zero expectation value.

Together with S6.156-S6.158, this completes the named
fixed-order Phi pole, b2, vacuum and stationary-H reference
assembly. It does not complete other coordinate dictionaries,
a global effective potential or physical truncation.
Finite bounds on J2 are not squared to infer higher-order
regulated vacuum-source terms.

## Independent verification

The report pins 18 source/proof/test files and 20 fields:
44 named identities, 44 scalar entries, 23 proof gates,
9 controls and 138 rejected inputs. Final private science
passed 214 tests in 22.23 seconds; fresh repository science
passed 214 in 22.58 seconds. Ordinary replay passed 239
tests in 2104.78 seconds. Independent native CLI replay
passed. The complete P8 snapshot passed 25306 tests in
3246.46 seconds with final exit code 0.

All 571 captured test files were present and unchanged.
Path-list SHA-256:

    dc0ec2bddbbcc7cdf440d8ef60ebeee9cee702e7be540b2959747a0ae9e985ed

The full-run adapter passed 128 original tuple comparisons.
Counters were 34095 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.160.
Only full regression uses the separately audited exact
GCD adapter; native, ordinary, CLI and direct science
retain unmodified SymPy with interpreter-only allowances.

The 34309-character native report was transferred losslessly
and all eighteen source hashes were verified. Report SHA-256:

    7bb4498c44c3a56b8e5b7b6a8c8c608193708c94bae80ce1b09e91db4f58a9a6

Independent tests compare fixed-counterterm vacuum
differentiation with the covariance source, reject moving
reference conditions, extract finite source products on
complex regulator circles and verify both light-line mass
derivatives and the endpoint estimates.

All cleanup preceded final private verification and freeze.
No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

Full derivative-action transport with the matched physical
source has passed native and fresh direct science and is
completing independent regressions. A complete second
low-energy elastic-cut allowance is in private development.
Neither is counted as published in this checkpoint.

Ordinary-Psi composite normalization is not inferred from
physical-source equivalence. Other required dictionaries,
physical finite-EFT truncation, V contours/cuts, finite-gravity
G and common-parent B remain open. No all-orders UV construction
is added to the adopted finite-EFT/necessary-positivity scope.
Scoped P8(a) and A.20-A.23 are unchanged; original P8(b) and
original P8 remain OPEN.
