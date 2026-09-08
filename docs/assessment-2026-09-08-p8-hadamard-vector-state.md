# P8: controlled Hadamard vector state

Original P8 remains open. The scoped P8(a) photon objective,
original linear classification and adopted V/G/B contract are
unchanged. This follows the
[clock-source audit](assessment-2026-09-08-p8-vector-clock-source.md).

## New result

[S6.55](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/hadamard/FORMULATION.md)
constructs a separately named all-order vector Cauchy state
without changing the action, physical clock or finite matching.
The frozen fourth-order state is not relabelled as Hadamard.

An exact all-order WKB recurrence and explicitly doubling,
time-independent momentum cutoffs define a locally finite
frequency and slope correction. A constrained Proca Hadamard
reference is built through an auxiliary globally hyperbolic
geometry that agrees with the physical one near the initial
Cauchy slice. An all-order comparison proves a smooth
two-point difference. The auxiliary metric is only a proof
device, not a replacement bounce.

The initial Bogoliubov change is bounded by 10^7/nu^6.
Exact physical readout identities transfer the energy,
pressure, their first derivatives and homogeneous clock-source
bounds to this new state. A separate complex-dimensional
dominated-convergence proof shows that this finite state
change adds no pole or finite counterterm.

At M*tau=10^12,m0*tau=1000, all five vector-only matched
bounds remain below 10^-14 of their reference scales.
This does not furnish other quantum fields, interacting loops,
complete functional response, corrected light cones, cutoff
or V/G/B closure.

## Verification

The report pins 15 sources and fully rebuilds S6.54 and its
ancestry. It checks 36 named identities comprising 61 scalar
entries, 26 proof checks and 113 invalid inputs.
The focused science suite passes 10 tests in 2.89 seconds;
the ordinary suite passes 46 tests in 347.19 seconds without
the broad GCD adapter. The separate read-only CLI passes.

A replay caught an invalid-input cache collision: SymPy
integers and floats could retrieve native-integer results
after warmup. Type-sensitive caches and a warm-cache regression
test corrected it before certificate generation.

The full P8 regression through S6.56 passes 5008 tests in
956.19 seconds, without excluding a frozen checkpoint. The
exact GCD adapter passes 128 original tuple comparisons and
records 5879 domain fallbacks and 6509 exact descents. This
run collected before the later S6.57 working tests were added.

The report SHA-256 is
8ba1aa16cee4dbe890b233fac864a83de569402e06b26c9c18ae56ef763a9200.
This is symbolic/rational verification and written continuous
and all-order proofs, not proof-assistant formalization.
