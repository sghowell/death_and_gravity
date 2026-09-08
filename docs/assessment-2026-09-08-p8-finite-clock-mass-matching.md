# P8: finite clock-mass vector matching and regulator limit

Original P8 remains open. The scoped P8(a) photon objective,
original linear classification and adopted physical-frame V/G/B
contract are unchanged. This follows the
[dimensional matching audit](assessment-2026-09-08-p8-dimensional-vector-matching.md).

## New result

[S6.53](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/clock_matching/FORMULATION.md)
specifies and verifies the scalar-coefficient continuation of the
actual clock-mass counterterms. It retains the entire frozen finite
potential, lapse jets, vector action and Gaussian preparation.
The new curvature terms are covariant on the original timelike
clock tube, linear in the vanishing mass deviations, and explicitly
use three times the spatial average in their dimensional extension.
This is a named first-variation prescription, not unique UV data.

Varying before clock restriction and before the regulator limit
reproduces all three actual lapse poles and the ordinary-Proca
finite control. The added mass terms have zero first homogeneous
spatial variation. At mu=m the matched second-order local energy
coefficient is

    -(10+alpha-3beta)H^2-(2/3)beta*H'.

The full fourth-order expression is retained exactly. Compact local
Euler variations of boxR vanish; no noncompact infinite-time flux
is erased and no finite scalar anomaly is discarded.

The continuation of the exact prepared mode pair is analytic on
the complex disk |D-3|<=1/4. Continuous rational coefficient bounds
control its exact/reference difference and full subtraction tail
uniformly over the time interval and all momenta. Their radial
tails are integrable uniformly in D. Dominated convergence gives
the physical finite integral at D=3; Morera's theorem gives
holomorphy. The pair does not conjugate complex D, and no
noninteger-dimensional Hilbert space is posited.

This supplies the previously missing limit needed to combine the
S6.51 finite mode integral with the matched local coefficients.
On u in [-1/2,1/2], at M*tau=10^12 and m0*tau=1000, each matched
vector energy/pressure magnitude is below 10^-14 of M^2/tau^2.
The error decomposition and rational continuous bounds are explicit.

## Remaining work

This is vector-only one-loop matching at a fixed background, not
a quantum-corrected solution. First time derivatives and the clock
source are under separate development in S6.54. All-order state
admissibility, higher functional variations, other field/higher
loops, corrected constraints/cones, interacting cutoff, a common
vacuum extension and finite-gravity V/G/B remain research requirements.
Small unsigned first-variation values do not protect a saturated
classical matter cone. No user-intervention blocker is present.

## Verification

The certificate pins 14 local sources and fully rebuilds S6.52 and
its frozen ancestry. It checks 43 named exact identities comprising
45 scalar entries, 24 proof checks and 59 rejected inputs.
Focused science passes **9 tests in 2.90 seconds**. The ordinary
suite passes **46 tests in 337.57 seconds**, without the broad GCD
adapter. The separate seeded read-only CLI passes.

The full P8 regression through S6.53 passes **4875 tests in 904.43
seconds**, without excluding a frozen checkpoint. The exact GCD
adapter passes 128 original tuple comparisons and records 5796
domain fallbacks and 6509 exact descents under the per-test seed
recipe. Source manifests and missing, extra or mutated report
fields are checked. The report SHA-256 is
17c6394918dd6bc1ca0238cab20f534db3950fb02a56eaac853e2d5377c9c1a4.

This is exact symbolic verification with written analytic proofs,
not proof-assistant formalization, peer review or original P8 closure.
