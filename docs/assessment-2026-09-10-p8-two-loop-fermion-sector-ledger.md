# P8: complete primitive two-loop fermion-sector ownership

Original P8 remains OPEN. This follows the
[spectral quadratic-sector audit](assessment-2026-09-10-p8-fermion-spectral-quadratic-sector.md).
No user intervention is needed for the current research.

## The complete primitive inventory

[S6.137](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/scheme/internal_momentum/spectral_bubble/light_pole/fermion_sector/FORMULATION.md)
integrates the quadratic fermions first, writing their
determinant as F=-Tr log D_F. The bosonic loop expansion
and Legendre transform give the genuinely mixed
order-two effective action

    Gamma_2,F=(1/2)Tr[(S_B^(0)'')^(-1) F''].

An explicit Gaussian calculation shows why the
connected two-fermion-loop dumbbell and the mixed
source/tadpole terms cancel in the Legendre transform.
The general regulated index argument is written out;
the finite-dimensional check is not used as a substitute
for that argument.

After exact Gaussian H elimination, the full scalar
Hessian is D^-1+W2, with W2 homogeneous quadratic in
Phi and retaining its nonlocal heavy kernel. Its
noncommuting inverse expansion gives six scalar rows
through four external Phi derivatives. The gauge
contraction supplies three further rows. In total
there are two vacuum, three two-point and four
four-point rows.

The color-singlet background makes the mixed Phi-A
fermion Hessian vanish. Literal SU(3) generators check
the open Casimir and all relevant traces. The active
scalar color/flavor factor is six, the active
contracted gauge factor eight, and the fourteen-flavor
gauge vacuum factor 56. Gauge Lorentz and vertex
factors remain intact; these color checks do not
evaluate a gauge integral.

## Counterterm ownership is separate from evaluation

The raw primitive trace is supplemented by bosonic
one-loop counterterm insertions, the fermionic
insertion -Tr[D_F^-1 D_F^(1)], and second-order
local references. Counterterms already paired inside
a renormalized proper subgraph cannot be counted
again in a separate list.

The ledger retains regulator-sensitive finite
products, including epsilon times a simple pole
and epsilon squared times a double pole. It also
retains the complete second-order field conversion
and the fundamental-cubic square G1^2+2G G2.

At this checkpoint's freeze, one paired quadratic
row was bounded, one quartic row had only its local
subset bounded, and seven primitive rows were
unevaluated. The tests reject falsely completed
or missing rows. A complete inventory is not a
complete error estimate.

## Independent verification

The native report pins 18 source/proof/test files
and 20 fields: 141 named identities and scalar
entries, 28 proof gates, 9 controls and 39 rejected
inputs.

Final private science: 268 tests in 21.06 seconds.
Independent repository science: 268 tests in 21.34 seconds.
Ordinary replay: 293 tests in 2358.05 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 20120 tests in 2880.18 seconds.

All 527 captured test files are present and unchanged.
Path-list SHA-256:

    1f8eb2d6578adb508ad31c70798aaf3b0a9126784f646fb5f1b58dd1a2dfe9ad

The adapter passed 128 original tuple comparisons.
Full-run counters: 34018 domain fallbacks, 7340
exact descents and 4 mixed fallbacks. The snapshot
predates S6.138's tests.

Only full regression uses the audited exact GCD
adapter. All native, ordinary, CLI and direct
science use unmodified SymPy with interpreter-only
allowances. The 35816-character native report was
transferred losslessly and all source hashes verified.
Report SHA-256:

    0f6a4e972f68a58287a6c3f31faf4b281bf57bc3ed7f0697f083a73141080f5e

Exact 22-file staging checks staged report/source
bytes and excludes nested continuations and unrelated
P4/P9 work. These are exact algebra and written
analytic proofs, not proof-assistant formalization
or independent peer review.

## Active continuation and nonclosure

S6.138 now bounds the full nonlocal heavy-vertex
quartic insertion family in an explicit local-parent
subtraction scheme. S6.139 supplies proper one-loop
fermion local references and finite field/coupling
conversion data. Both are native-certified with
fresh replays underway. Neither is being promoted
to the complete common two-loop matching result.

Remaining primitive integrals, proper-counterterm
insertions, finite interaction conversion and
higher-loop/truncation control remain research
obligations. The adopted V contours, finite-gravity
G and common bounce-parent B are still open.
No all-orders UV construction is added to the
agreed finite-EFT/necessary-positivity task.
Scoped P8(a) and A.20-A.23 are unchanged. Original
P8(b) and P8 remain OPEN.
