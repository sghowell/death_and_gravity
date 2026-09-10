# P8: complete one-loop vacuum forward coefficient

Original P8 remains OPEN. This follows the
[one-loop light-pole audit](assessment-2026-09-09-p8-one-loop-light-pole.md).
No user intervention is needed for the current research.

## Completed full one-loop bound

[S6.113](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/FORMULATION.md)
bounds the actual polynomial vacuum's complete one-loop forward
four-point coefficient, with the mass-one light pole and unit residue
fixed by S6.112. The exact heavy kernel is retained throughout loop
momentum integration. The calculation includes the mixed terms and
all three scattering channels.

A literal finite-site quartic Hessian and five further derivative
checks fix the full nonlocal vertex and bubble symmetry factor.
Complex on-shell external momenta, their permissible contour shift
and a positive Feynman-parameter denominator establish the analytic
disc used to extract the coefficient. No expansion in loop momentum
divided by the heavy mass is made.

Writing M for heavy mass squared, g for cubic coupling squared and
c=g/M^2, the exact bound is

    |delta b2 at one loop| < 60 c^2 M + 14000 c^2
                          < 10^-6 times 4 lambda.

The first term bounds the complete angular/shift remainder; the
second bounds the twice-differentiated ultraviolet-subtracted radial
part. The t-channel remainder is included. All local counterterms,
the once-fixed finite potential contact and external normalization
are accounted for without introducing a free b2 subtraction.

Since the tree coefficient is exactly 4 lambda, the tree-plus-one-loop
coefficient is strictly positive. This is a finite-order result,
not an all-order positivity theorem or an all-energy UV completion.

## Independent verification

The lossless native report pins 19 source/proof/test files and all
22 report fields. It contains 54 named exact identities, 69 scalar
entries, 44 proof gates and 16 rejected inputs.

Final private science: 140 tests in 10.12 seconds.
Independent repository science: 140 tests in 10.06 seconds.
Ordinary read-only replay: 167 tests in 2486.63 seconds.
Independent native CLI: passed.
Complete P8 snapshot: 10778 tests in 3355.66 seconds.

All 479 captured test files are present and unchanged.
Path-list SHA-256:

    9278c31f0ae09b52ee80b6fc50821656e0a9756439dedd81b01dd41868af3d5e

The unchanged exact GCD runner passes 128 original tuple comparisons.
Full-run counters: 32785 domain fallbacks, 7340 exact descents and
4 mixed fallbacks. The snapshot predates S6.114's two test files.
Native generation, ordinary replay and CLI use unmodified scientific
SymPy; the separately audited adapter is confined to full regression.

Native generation reused the S6.112 batch's byte-identical cached
ancestry only after verifying every previously captured P8 Python,
Markdown and JSON input hash. The child's ordinary, CLI and full
snapshot replay each ran independently in a fresh process.

Report SHA-256:

    336e83ee3053018d9ce8a84a90651e0ceab79eb2d176003c5e466faf126fc2c2

Exact 23-file staging verifies the report and all 19 staged source
hashes, plus only the three named root audit files. Nested children
and unrelated P4/P9 changes are excluded. Analytic contour,
integration and perturbative estimates are source-pinned written
proofs, not proof-assistant formalized or peer reviewed.

## Continuation and open obligations

The elastic-cut child has passed its science, ordinary, CLI and full
regressions and awaits its own publication. The spin-two form-factor
child has a native report and is in independent replay. Spectral-cut,
reference-flow and full quantum field-map work continue.

The current positive coefficient does not supply higher-loop
remainders, high-energy contour boundedness or a finite-gravity
IR/Regge subtraction. Nor does it identify this flat-space polynomial
vacuum with a controlled propagating cosmological bounce parent.
Those V/G/B requirements, absolute coupled renormalization and
corrected cosmological cones remain open. Old states and counterterms
are not transferred.

The scoped P8(a) objective and A.20-A.23 extensions remain complete.
Original P8(b), and therefore original P8, is not finished or closed.
