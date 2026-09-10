# P8: complete first fermion four-scalar increment

Original P8 remains OPEN. This follows the
[local-reference matching audit](assessment-2026-09-10-p8-fermion-local-reference-matching.md).
No user intervention is needed for the current research.

## Complete finite-mass box, not only a local threshold

[S6.132](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/FORMULATION.md)
computes the complete first fermion contribution to the four-scalar
1PI function at the named GY14 boundary. Both Yukawa signs, all active
colors/flavors, six cyclic box orderings and both orientations are
retained. The inert flavors have no Phi four-point vertex.

An independent Dirac trace fixes the zero-momentum finite vertex
v4=-64NY^2/Q and its whole ultraviolet pole 24NY^2 Ibar/Q,
with N=6 and Q=16pi^2. Scattering has the opposite sign to
the Euclidean 1PI vertex.

The constant-background kinetic function fixes the degree-two vertex.
A complete permutation-invariance calculation on the ten quadratic
momentum contractions has rank eight; momentum conservation reduces
the two orbit sums to the single allowed quadratic operator.
Its on-shell value is constant, so it has no forward b2.

The proof does not stop at this local expansion. The six boxes'
entire higher-momentum sum is bounded on the complex forward disc
of radius five. Separate evenness in energy and spatial momentum
removes square-root coordinate artifacts. Dirac resolvent norms,
radial moments and a convergent full composition sum give

    E_box < 3 * 10^8 Y^2/(Q mF^4).

This bounds every momentum degree at this one-loop order.
Cauchy's estimate controls the forward second coefficient; it is
not a claim to sum every loop order or justify a high-energy contour.

## Canonical references and no double counting

With the S6.131 fermion kinetic factor kappa=1-f'(1), the selected
scalar-tree plus fermion amplitude is (A0+A_F)/kappa^2.
If the finite potential quartic shift v4 is placed in the tree
potential, the residual box must instead be A_F+v4.
The report verifies this exact dictionary and its first-order form.

At the actual boundary the box-error ratio to the tree b2 is below
10^-604. The fermion residue contribution is retained with its sign.
The explicitly selected tree-plus-fermion normalized interval is
positive; a generic valid but wider error interval can be inconclusive.
This result does not yet include the old scalar loop's finite scheme
conversion. S6.133 supplies that separate complete one-loop combination.

## Independent verification

The native report pins 18 source/proof/test files and 20 fields.
It contains 127 named identities, 254 scalar entries, 29 proof gates,
9 controls and 43 rejected inputs.

Final private science: 247 tests in 20.34 seconds.
Independent repository science: 247 tests in 20.55 seconds.
Ordinary replay: 272 tests in 2112.68 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 18672 tests in 2895.83 seconds.

All 517 captured test files are present and unchanged. Path-list SHA-256:

    938164bb53f1fa6d57d34fa8e488049181e048f468ae9068be85fb3c6b785d3f

The adapter passed 128 original tuple comparisons.
Full-run counters: 33983 domain fallbacks, 7340 exact descents,
4 mixed fallbacks. This captured snapshot predates S6.133's tests.
Only the full regression uses the audited exact GCD adapter; native,
ordinary, CLI and both direct science runs use unmodified SymPy.

The first native export hit Python's 4300-digit integer-string limit.
It did not produce a saved report and is not counted as successful.
The scientific files were unchanged; a fresh native ancestry rebuild
used the documented interpreter-only allowances (recursion 4000 and
unlimited formatting of these trusted exact integers). The
[durable replay wrapper](p8-replay-runtime.md) reproduces those settings.

The 97514-character native report was transferred losslessly and every
source hash verified. Report SHA-256:

    abbc5dfd0adffad21deb187280c1da862845c838ab023c87c8713166320f935d

Exact 22-file staging contains this checkpoint, this root audit,
CLAIMS.md and README.md only, and verifies staged source/report bytes.
It excludes nested continuation checkpoints and unrelated P4/P9 work.
The analytic arguments are written proofs, not proof-assistant
formalization or independent peer review.

## Active continuation and original obligations

S6.133 combines the scalar and fermion one-loop contributions in
one explicit canonical reference scheme; its direct science,
ordinary and CLI replays passed. S6.134 extends the insertion bound
to unbounded internal momentum and an explicitly finite reference
energy window; its native report and both direct science runs passed,
with fresh replays running.

The next integration uses a directly derived positive spectral
representation of the fermion insertion, with inner reference
pairing and an overall local subtraction. Each new result is
limited to its proved contribution family.

New-model later-loop errors, justified V contours, finite-gravity
IR/Regge allowance, and common-parent bounce field/state/cutoff
control remain unproved. These are the adopted finite-EFT and
necessary-positivity obligations, not a new demand to construct an
all-orders UV theory. Scoped P8(a) and A.20-A.23 remain unchanged.
Original P8(b) and P8 remain open.
