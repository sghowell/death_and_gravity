# P8: quantitative fifth-order vector-stress regularity

Original P8 remains open. This follows the
[actual vector variation audit](assessment-2026-09-08-p8-actual-vector-variation.md).
The named action, exact vector state, fixed Cauchy cutoffs and
finite matching prescription are unchanged.

## New result

[S6.59](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/regularity/FORMULATION.md)
proves explicit continuous bounds for the finite vector energy
and pressure through five physical-time derivatives on
u in [-1/2,1/2], for every comoving momentum and m0*tau>=1000.

An eighth-order WKB reference is used only for comparison.
Its actual residual and derivative are bounded exactly for
both polarizations. The unchanged state preparation is controlled
in three frequency bands determined by its existing cutoffs.
The oscillatory evolution estimate retains the nonzero initial
mixing and both integration-by-parts endpoints.

The physical differentiated readouts are derived with the complete
fixed-comoving chain rule. Exact-ODE rows evaluated on reference
products are distinguished from derivatives of an approximate mode.
For all four readouts and all six derivative orders, the necessary
low subtraction coefficients vanish exactly. A lower-reference
failure control detects a nonzero coefficient at second order;
it is not misreported as divergence of the exact state's derivative.

Continuous radial bounds then control the full momentum integral.
The actual finite local matching terms are differentiated as well.
At M*tau=10^24,m0*tau=1000, all twelve normalized bounds are
below 10^-18. The fifth-order energy and pressure bounds are
below 1.165*10^-19 and 1.132*10^-19, respectively. The j-th
physical derivative is normalized by M^2/tau^(2+j).
This additional scale example does not replace the earlier
M*tau=10^12 results.

## Remaining work

These are fixed-background C5 estimates, not metric-functional
response bounds or a self-consistent quantum solution. The
contact/retarded second variation, compatible initial covariance
response, corrected coupled cones, interactions, cutoff and V/G/B
remain open.

An additional off-clock vector-block check motivates a separately
named, selected-state scalar-profile candidate, now being tested.
Its purpose is to preserve the original clock history in the
retained Gaussian sector with explicit coefficient bounds. No
stability, full quantum or UV conclusion is inferred from that
background-level construction. There is no user-intervention blocker.

## Verification

The report pins 14 sources and fully rebuilds S6.58. It verifies
86 named identities comprising 381 scalar entries, 98 continuous
proof checks and 120 rejected inputs. Focused science passes
14 tests in 140.37 seconds. The ordinary suite passes 50 tests
in 525.88 seconds without the broad GCD adapter; the independent
read-only CLI passes.

The full regression through S6.59 passes 5148 tests in
1048.64 seconds without excluding a frozen checkpoint. The
adapter passes 128 original tuple comparisons and records
6458 domain fallbacks and 6509 exact descents. This run
collected before the S6.60 working tests were added.

Report SHA-256:
ee3da10ed767ca94be285bf262b745d12f8a3b0c5dfa3b2a7797da433c2c58a2.
The evidence is exact symbolic/rational computation with written
continuous proofs, not proof-assistant formalization.
