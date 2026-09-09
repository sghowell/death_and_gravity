# P8: new ordinary-Proca finite local response

Original P8 remains open. The completed scoped P8(a) objective is unchanged,
and no user intervention is needed to continue. This follows the
[central scalar-cone audit](assessment-2026-09-08-p8-central-scalar-cone.md).

## New result

[S6.84](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/FORMULATION.md) rebuilds the finite local
homogeneous response for the source-free constant-physical-mass Proca
candidate. The old mass-source response is not transferred.

All mass jets are set to zero before the actual clock substitution.
The common dimensional calculation, including the derivative of its
counterterm before D=3, agrees independently with the covariant finite
heat action. In the normalization 1/(64*pi^2), its density is

    (5/2)m^4+(5/3)m^2 R-4 a4_scalar.

Both coordinate Euler currents and normalized physical energy/pressure
variations are calculated. Constant spatial-rescaling terms cancel
exactly in the physical stress. Weighted self-adjointness is checked
for coordinate Hessians, not incorrectly imposed on normalized stress.

The actual nonlinear clock-metric pullback keeps its second-map contact.
Only the already-fixed local part of the S6.82 profile is included here.
Its order-zero potential cancels exactly. No state, profile or finite
counterterm is selected again in response to a source.

The finite fourth-order clock-chart coefficient is -4 e e^T with
e=(1/(2h),1). It has rank one and null vector (1,-1/(2h)).
This is not the old invertible two-source block and is not, by itself,
a statement about a resummed spectrum or full quantum characteristics.

Continuous coefficient enclosures on I=[-1/2,1/2] give local C4-to-C0
response bounds below 10^-790 at fixed m=1000,L=10^400, in both
physical and clock coordinates. Nonlocal response, the full profile,
no-loss inversion and arbitrary spatial sources are outside this
local checkpoint.

## Verification

The losslessly saved native report pins 16 source files and fully rebuilds
S6.83. It checks 246 named identities, 250 scalar entries, 15 proof gates
and 100 rejected inputs.

The native science suite passes 30 tests in 23.50 seconds. The ordinary
native suite passes 60 tests in 1963.57 seconds, and the independent
read-only CLI passes. These use unchanged native scientific arithmetic;
focused pytest runs use the audited static namespace collection setup.

The complete P8 regression passes 6452 tests in 2305.90 seconds.
All 413 captured files are present and unchanged, with path-list SHA-256

    af05bfd5f47483d8a73315dc8d5101cdeebe2c4f720abbbc3f27d476f2b19884

The unchanged exact GCD adapter passes 128 original tuple comparisons
and records 27297 domain fallbacks, 7340 exact descents and four mixed
fallbacks. The snapshot predates S6.85's added tests.

Report SHA-256:

    240b1a0550969a51bd1c0ea4043360566a1b1cd10408ea30b2f89d560dc1aea2

This is exact algebra and continuous-domain analysis with written proofs,
not proof-assistant formalization.

The final whitespace check reports only extra blank EOF lines in the seven
frozen formulation/proof Markdown files and four frozen verifier/test
files. Their certified bytes are preserved. Ruff passes and no other
whitespace warning is present.

## Continuing work

S6.85's native report has been saved losslessly with all 19 source hashes
verified; its 30 independent science tests pass and the ordinary, CLI and
full replays are running. The new isolated scalar inverse has passed its
30 native science tests and its sources are frozen for report generation.
The actual coupled-inverse draft has passed 37 independent temporary
tests and is awaiting its parent/report pipeline.

Quantum stability/cones, nonlinear and spatial control, canonical profile
interactions, omitted operators, mixed loops, thresholds, finite Wilson
matching and common-parent V/G/B requirements remain open. No UV
completion, universal exclusion or original P8 closure is asserted.
