# P8: actual one-loop low-energy elastic cut

Original P8 remains OPEN. This follows the
[full forward-coefficient audit](assessment-2026-09-09-p8-one-loop-forward-coefficient.md).
No user intervention is needed for the current research.

## Completed elastic-cut calculation

[S6.114](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/FORMULATION.md)
computes the actual polynomial vacuum's one-loop two-light-particle
discontinuity for mass-one invariant s between four and six. It uses
the complete heavy-exchange tree vertex, not its derivative expansion.

An anchored exact angular primitive retains contact/exchange
interference and has a checked continuous threshold limit. Independent
phase-space and constant-bubble calculations fix both the identical
intermediate-particle factor and the factor relating 2 Im A to Im A.
Physical-angle and invariant monotonicity give bounds over the entire
window, including a strictly positive subwindow.

With rho(s)=Im A at one loop, the explicitly computed partial
forward-cut contribution satisfies

    lambda^2/20 < (2/pi) integral_4^6 rho(s)/(s-2)^3 ds
                < 3 lambda^2.

Subtracting this cut contribution and S6.113's complete one-loop
real-coefficient error from the tree value 4 lambda still leaves a
positive coefficient, with combined loss below one millionth of
that tree value.

This is a finite-order partial improved functional. A global
dispersion identity, control of its high-energy contour, higher-loop
remainders and a uniform heavy-resonance treatment are not assumed.
The exact threshold density is zero; strict positivity is asserted
above threshold, not at it.

## Independent verification

The lossless native report pins 17 source/proof/test files and all
21 report fields. It contains 33 named exact identities and scalar
entries, 32 proof gates and 16 rejected inputs.

Final private science: 107 tests in 9.82 seconds.
Independent repository science: 107 tests in 11.60 seconds.
Ordinary read-only replay: 133 tests in 2442.83 seconds.
Independent native CLI: passed.
Complete P8 snapshot: 10911 tests in 3347.73 seconds.

All 481 captured test files are present and unchanged.
Path-list SHA-256:

    34167e45fe29f2cdbd422bc4a28c2a1971e805e610bcff86c46cc224e3ca760d

The unchanged exact GCD runner passes 128 original tuple comparisons.
Full-run counters: 32770 domain fallbacks, 7340 exact descents and
4 mixed fallbacks. The snapshot predates S6.115's two test files.
Native generation, ordinary replay and CLI use unmodified scientific
SymPy; the separately audited adapter is confined to full regression.

Native generation reused the same batch's byte-identical cached
ancestry only after checking every previously captured P8 Python,
Markdown and JSON input hash. This child's ordinary, CLI and full
snapshot replay each ran independently in a fresh process.

Report SHA-256:

    558e10ebc33b3a4ff9a0b09a9acf479cd18e1d6f2c2b62d014c0c0b855701f12

Exact 21-file staging verifies the report and all 17 staged source
hashes, plus only the three named root audit files. Nested children
and unrelated P4/P9 changes are excluded. Continuous cut inequalities
and perturbative unitarity arguments are source-pinned written proofs,
not proof-assistant formalized or peer reviewed.

## Active continuation

The spin-two form-factor child has passed private/repository science
and native report generation; independent replays are running. It
computes only the t-channel gravitational matter-vertex correction,
not the full finite-gravity remainder. Its spectral-cut successor is
rebuilding native ancestry. The reference-flow successor has passed
private science. Work on the full quantum field map now includes a
literal derivative-polynomial Wick check with its finite evanescent
subtraction term retained.

Higher-loop and all-energy vacuum control, finite-gravity IR/Regge
bounds, absolute coupled renormalization, corrected cosmological
cones and a controlled common propagating bounce parent remain open.
No old state or counterterm is transferred.

The scoped P8(a) objective and A.20-A.23 extensions remain complete.
Original P8(b), and therefore original P8, is not finished or closed.
