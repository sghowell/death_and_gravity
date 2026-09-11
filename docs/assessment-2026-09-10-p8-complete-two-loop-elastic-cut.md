# P8: complete two-loop low-energy Phi elastic cut

Original P8 remains OPEN. This follows the
[physical-source transport checkpoint](assessment-2026-09-10-p8-two-loop-physical-source-transport.md).
No user intervention is required for the current research.

## The complete first amplitude on the physical cut window

[S6.161](../problems/P8/s6/continuation/s6_161/FORMULATION.md)
bounds the complete canonical first-loop GY14 Phi amplitude
on physical s in [4,6], not just its forward Taylor coefficient.
It retains all three scalar channels, full heavy kernels,
MS scale/contact terms, the full fermion local reference
and convergent Dirac difference, and the sole first
canonical field correction -2 k0 A0.

The scalar bubble logarithm is integrated across its
continuous two-light threshold. Triangle and box parameter
integrals are treated as Feynman boundary values: change
variable through the simple zero, subtract the smooth
weight at the zero, and only then bound the result.
For the box, integration by parts at nonzero regulator
reduces the double-pole expression before this subtraction.
No absolute integral through an unsplit double pole is used.

For heavy mass squared M>=32 the resulting parameter bounds are

    triangle <=(3 log M+22)/M,
    box      <=(6 log M+52)/M^2.

A proved rational logarithm cap is used at the actual M.
The fermion estimate retains the entire zero-momentum MS
reference 64 N Y^2/Q and bounds the remaining momentum
dependence using all six cyclic Dirac words. It is a bound
on the full amplitude, not an extrapolation of a b2 bound.

The actual complete first-amplitude upper bound is

    B1=5.4147127665826475859e-408 <1e-405.

It need not be small relative to the exceptionally canceled
tree amplitude. The cut calculation uses its absolute norm.

## Complete second cut and improved formal margin

With physical Phi mass fixed at one, only the two-Phi
intermediate-state tree/one-loop interference contributes
through loop two on this window. Three-Phi states are
forbidden by parity. Heavy and fermion thresholds lie above
the window; gauge and four-Phi contributions first enter
at later loop order. No gauge confinement gap is assumed.

Writing beta=sqrt(1-4/s), the identical-particle normalization is

    rho1=beta/(32pi) integral_0^1 A0^2 dz,
    rho2=beta/(16pi) integral_0^1 A0 Re(A1) dz.

The factor two in the interference is retained. The known
tree bound |A0|<73 lambda and the dispersion weight on [4,6]
give

    |I2| <=73 lambda B1/1280,
    |I2|/(4lambda) <1e-407.

Numerically the latter is
7.7201959367291655033e-410. The earlier
lambda^2/20<I1<3lambda^2 bound is unchanged.
Subtracting both computed cut orders from the complete
through-two-loop b2 polynomial preserves a positive
formal uniform lower bound; the total relative allowance
remains below 1e-6.

No sign is imposed on rho2. The square of A1 is only one
piece of the later three-loop cut and is not mislabeled as
that complete cut. A positive truncated margin does not
bound the physical omitted remainder or establish an exact
global dispersive equality.

## Independent verification

The native report pins 18 source/proof/test files and 20
fields: 52 named identities, 52 scalar entries, 23 proof
gates, 9 controls and 234 rejected inputs. Final private
science passed 333 tests in 22.06 seconds; fresh repository
science passed 333 in 21.99 seconds. Ordinary replay passed
358 tests in 2170.65 seconds. Independent native CLI replay
passed. The complete captured P8 snapshot passed 25910
tests in 3253.82 seconds with final exit code 0.

All 575 captured test files were present and unchanged.
Path-list SHA-256:

    804fd5f13b8d337434223df2688b9e97ed192292179d522716f8942a7d1d75ac

The full-run adapter passed 128 original tuple comparisons.
Counters were 34105 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.162.
Only full regression uses the separately audited exact
GCD adapter; native, direct science, ordinary and CLI retain
unmodified SymPy with interpreter-only allowances.

The 71155-character native report was transferred losslessly
in six chunks; all eighteen source hashes were independently
verified. Report SHA-256:

    2ef5a0230f712fc5f954ed970e67cef74b86b4216be659075a07b414357d5f13

Independent tests compare physical threshold logarithms,
boundary-subtracted triangles/boxes, analytic radial box
forms, direct Euclidean parameter integrals, explicit Dirac
matrix singular values and optical-theorem factors.
The written proof records its primary-source checks and
separates physical positive-state counting from a
finite-order amplitude calculation.

Cleanup preceded final private verification and byte freeze.
No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

S6.162 has passed native and fresh direct science for the
full finite-kappa analytic target's restricted classical
match. Its independent regressions are running; it is not
counted as published here. A separately named saturated
Yukawa profile is being tested privately for pointwise
mass protection and preservation of named low-order data.

Quantum target matching, physical finite-EFT truncation,
global V contours/cuts, finite-gravity G and common-parent
B remain open. A mass floor alone is not a rolling-state
or derivative-map dictionary. The adopted finite-EFT and
necessary-positivity scope is unchanged. Scoped P8(a) and
A.20-A.23 are unchanged; original P8(b) and original P8
remain OPEN.
