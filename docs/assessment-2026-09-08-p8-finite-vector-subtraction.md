# P8: a finite prescribed vector energy and pressure subtraction

Original P8 remains open. The scoped P8(a) photon objective, original
linear CD/M1 classification and physical-frame V/G/B contract are
unchanged. This follows the [physical energy comparison](assessment-2026-09-08-p8-vector-energy-comparison.md).

## New result

[S6.51](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/subtraction/FORMULATION.md)
supplies a complete order-zero, two and four adiabatic subtraction for
the actual retained vector's physical energy and pressure. It retains
the clock mass jets, exact canonical normalization and the selected
S6.50 Gaussian preparation on the actual rolling interval. A separate
formal derivative parameter determines the subtraction order; finite
parts of those complete coefficients are not silently discarded.

The reference remainder is an exact rational identity, not an
unbounded asymptotic series. The continuous coefficient bounds give
at most D/(4*a_s^3*omega^5) per polarization, D=400000, for both
the absolute energy and pressure remainders. The full R^3 momentum
integral, with two transverse and one longitudinal mode, is at most
D/(72*m^2), m=m0*tau>=1000. No momentum cutoff is needed.

Combining this with the independently bounded exact-mode/reference
difference gives a finite, explicitly prescribed subtraction for the
exact Gaussian observables, uniformly on u in [-1/2,1/2]:

    abs(rho_sub)<=12C/m+D/(72m^2),
    abs(p_sub)<=108C/(5m)+D/(72m^2), C=2000000.

Physical densities are tau^-4 times these normalized expressions.
Relative to M^2/tau^2, at M*tau=10^12,m0*tau=1000 the bounds
are 4320001/(180*10^24) and 7776001/(180*10^24). The reference
tail alone contributes at most 1/(180*10^24). Both integrand
differences in this decomposition are absolutely integrable.

An independent ordinary-Proca transverse calculation reproduces the
integrated second- and fourth-order benchmark. For the actual
clock-sensitive energy, stopping at second order leaves a nonzero
logarithmic UV divergence at the bounce: the fourth-order coefficient
tends to -4/27, giving -log(Lambda)/(54*pi^2). The finite exact
evolution difference cannot cancel it.

## Boundary and next work

The result is a complete finite integral in this stated subtraction
prescription, not a scheme-independent bound or a full covariantly
matched quantum correction. The finite mass/metric counterterms and
all required state admissibility still need matching. In particular
finite-order initial data are not called all-order Hadamard. Other
field loops, all-time state tails, corrected constraints/cones, cutoff
and V/G/B are not supplied. No frozen action or earlier matched
finite counterterm has changed; no user intervention is needed now.

The next calculation retains dimensional polarization counts and
canonical friction through integration. It checks local poles and
finite terms against counterterm variations. Its ordinary-Proca
diagnostic already shows why taking a component-wise finite part
without first varying the dimensional counterterm action is insufficient.
The actual clock-dependent finite extension must preserve S6.47's
already frozen finite potential and lapse jets.

## Verification

The report pins 10 local sources and fully rebuilds S6.50 and its
ancestry. It checks 11 named exact scalar identities, 13 continuous
proof checks and 50 rejected inputs. The scientific suite passes
**11 tests in 83.25 seconds**. The ordinary suite passes **42 tests
in 430.61 seconds**, without the broad GCD adapter; the separate
seeded read-only CLI also passes.

The full current P8 regression passes **4792 tests in 879.08 seconds**,
including every frozen checkpoint and the four then-collected working
dimensional science tests. The exact GCD adapter passes 128 original
tuple comparisons and records 5763 domain fallbacks and 6509 exact
descents under the per-test seed recipe. No frozen checkpoint is
excluded. Source manifests and missing, extra or mutated report fields
are checked. This is exact symbolic verification plus written continuous
proofs, not proof-assistant formalization or external peer review.
Original P8 is not finished or closed.
