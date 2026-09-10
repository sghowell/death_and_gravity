# P8: proper fermion local-reference matching

Original P8 remains OPEN. This follows the
[full heavy-vertex insertion family](assessment-2026-09-10-p8-full-heavy-vertex-fermion-insertion.md).
No user intervention is needed for the current research.

## Proper fermion references at the shared MS scale

[S6.139](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/scheme/internal_momentum/spectral_bubble/light_pole/fermion_sector/heavy_vertices/fermion_references/FORMULATION.md)
derives the one-loop proper fermion kinetic,
mass and zero-momentum Phi-vertex references
for scalar and Feynman-gauge exchange.
The Euclidean inverse-resolvent derivation
fixes the proper self-energy signs. Literal
Clifford and color algebra retains the
open-line Casimir Cf=4/3, distinct from a
closed fermion-loop multiplicity.

Dimensional numerator terms are retained
until after MS pole subtraction. The gauge
mass and kinetic kernels contain the finite
remainders -2 and -1 respectively. At zero
Euclidean momentum their normalized finite
values are 2 and 1/2.

For scalar mass ratio r in [0,1], the exact
parameter integrals J0,J1,R have continuous
massless and equal-mass endpoint values.
The proper Phi vertex is the mass derivative
at fixed renormalization scale; the scale is
set to m only afterwards. Resetting it before
differentiation would give the wrong gauge
finite vertex. The exact shared coefficients are

    z=(-Y J1+a Cf/2)/Q
    eta=(Y J0+2a Cf)/Q
    ups=[Y(J0+2R)-6a Cf]/Q.

The ultraviolet counterterms independently
reproduce the frozen one-loop beta function
y(15Y-8a)/Q. The twelve inactive Yukawa flavors
have gauge kinetic/mass corrections but exactly
zero Phi Yukawa vertex.

## Finite maps and explicit remaining work

The same finite Phi normalization w from
S6.133 is retained. At the selected order,

    m_zero/m_MS=(1+h eta)/(1+h z)
    y_zero/y_MS=(1+h ups)/
       [(1+h z) sqrt(1+h w)].

These are reference maps, not an assertion
that MS input parameters remain unchanged
under canonical normalization. The selected
mass and Yukawa relative bounds are about
4.31508193781e-207 and 2.79862203716e-206.
The code rejects inappropriate exact domains
and does not claim positivity for large
uncontrolled input corrections.

This supplies local proper counterterms and
finite anchors for subsequent insertions.
It does not evaluate those insertion integrals
or the remaining primitive rows. Finite
S6.138 outer-reference conversion remains
separate, as do old-box field/parameter
conversions and higher regulator orders when
later forest products require them.

## Independent verification

The native report pins 18 source/proof/test
files and 20 fields: 76 named identities,
357 scalar entries, 30 proof gates, 9 controls
and 68 rejected inputs.

Final private science: 209 tests in 20.55 seconds.
Independent repository science: 209 tests in 20.51 seconds.
Ordinary replay: 234 tests in 2362.53 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 20646 tests in 2895.71 seconds.

All 531 captured test files are present and unchanged.
Path-list SHA-256:

    d1bcd332e5a2ca39b004695aa75e56959d70fce5f4c24967b169d7df88d22357

The adapter passed 128 original tuple comparisons.
Full-run counters: 34036 domain fallbacks,
7340 exact descents and 4 mixed fallbacks.
The snapshot predates S6.140's tests.

Only full regression uses the audited exact GCD
adapter. Native, ordinary, CLI and direct science
use unmodified SymPy with interpreter-only
allowances. The 48696-character native report
was transferred losslessly and all source hashes
verified. Report SHA-256:

    f7ae3eeb624a0d733eb0a0d27fd7b17c56aec0e393c3cebad9021a0947f351a8

Exact 22-file staging checks staged source/report
bytes and excludes nested continuations and
unrelated P4/P9 changes. The analytic arguments
are written proofs, not proof-assistant
formalization or independent peer review.

## Current continuation and nonclosure

S6.140 is native-certified, with fresh ordinary,
CLI and full replays running. It bounds the
24 self-energy-chord words per scalar/gauge
quartic primitive. S6.141's 12 opposite-chord
words have passing private and repository
science; a fresh native replay is underway.

S6.141 now lives in a shorter continuation
directory after macOS rejected an absolute
test path of length 1025. All 18 scientific
files were relocated byte-for-byte unchanged;
the failed-path copy is retained in a temporary
archive. No older scientific source was moved
or edited. The source content, not the failed
import, is what the new replay verifies.

The remaining 24 vertex words per sector are
under active research. Other primitive rows,
forests, finite interaction/field conversion,
full two-loop error and truncation control
remain open. V contours, finite-gravity G and
common bounce-parent B remain open. No
all-orders UV construction is added to the
adopted finite-EFT/necessary-positivity task.
Scoped P8(a) and A.20-A.23 are unchanged.
Original P8(b) and P8 remain OPEN.
