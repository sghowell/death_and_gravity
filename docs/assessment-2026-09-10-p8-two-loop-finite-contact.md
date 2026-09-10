# P8: once-fixed finite potential contact at two loops

Original P8 remains OPEN. This follows the
[wineglass audit](assessment-2026-09-10-p8-two-loop-wineglass.md).
No user intervention is needed for the current research.

## Completed fixed-condition contact insertion

[S6.125](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/FORMULATION.md)
evaluates the finite quartic contact already fixed by the
constant-field potential condition in S6.110 and retained
in S6.113. It is not a new parameter fitted to the
forward scattering coefficient.

The literal two-field Hessian is differentiated before
restricting to the classical stationary heavy field.
Its Schur kernel F is strictly between zero and its
high-momentum limit A. After the entire I0 reference
subtraction, the remaining finite contact is strictly
negative. Two independently anchored heavy-mass radial
integrals determine its exact value. Its magnitude is
about 1.0046945312 times 10^-409 and is bounded by
(50/3) lambda4^2.

Although this contact has zero tree-level second momentum
derivative, its insertion in the one-loop amplitude does
not vanish. At fixed g and M it is the quartic variation
of the complete renormalized one-loop amplitude. Both
heavy triangles and the linear variation of the inherited
local counterterms remain present.

All channel weights and the integrated variation give

    E_contact = (10025/18) lambda4^3 < 7 times 10^-612,
    E_contact / (4 lambda) < 2 times 10^-12.

Adding the four disjoint raw-graph groups gives known
two-loop contributions below 3 times 10^-607 and 10^-7
of tree. The known one-plus-two-loop contributions
remain below 10^-6 of tree. Full two-loop normalization
and source-aware matching are not inferred from this sum.

## Independent verification

The native report pins 18 source/proof/test files and 20 fields:
40 named exact identities and scalar entries, 34 proof gates,
9 controls and 14 rejected inputs.

Final private science: 111 tests in 24.32 seconds.
Independent repository science: 111 tests in 24.27 seconds.
Ordinary read-only replay: 136 tests in 2178.20 seconds.
Independent native CLI: passed.
Complete P8 snapshot: 15311 tests in 2966.50 seconds.

All 503 captured test files are present and unchanged.
Path-list SHA-256:

    3b977930748a39de36711edad9b326ef81dc153d24e7319e21081f56aabbb914

The audited exact GCD adapter passes 128 original tuple comparisons.
Full-run counters: 33969 domain fallbacks, 7340 exact descents,
4 mixed fallbacks. The snapshot predates S6.126's two test files.
Native generation, ordinary replay and CLI use unmodified
scientific SymPy; the adapter is confined to full regression.
Ordinary, CLI and full regression each ran in a fresh process.
Native ancestry reuse followed immutable P8 input-hash checks.

The 54772-character native report was transferred losslessly
and all 18 source hashes were verified. Report SHA-256:

    80ea67b9827236b172bc7f639988d293f64eac9a7b27de0d80533f7afcb4ace8

Exact 22-file staging checks every staged source/report hash
and includes only this checkpoint, this root audit, CLAIMS.md
and README.md. Nested children and unrelated P4/P9 work are excluded.

The continuous sign, integral, local subtraction and Cauchy
arguments are source-pinned written proofs, not proof-assistant
formalized or peer reviewed. A small finite contact is not
a statement about bare ultraviolet stability.

## Active continuation and remaining scope

The canonical two-loop light-pole child is now frozen and
native-certified: 743 exact identities, 32 gates, 33 rejected
inputs, and 926 private plus 926 repository science tests.
Its report and all 21 source hashes were verified.
Fresh ordinary, CLI and complete regressions are running.
It derives a unit-disc quadratic self-energy upper bound
below 10^-18, with the fixed-order mass-one pole and unit
residue. This is not an all-orders propagator statement.

The next analysis checks the independent complete canonical
two-loop four-point counterterm and external-field ledger.
A two-loop source-aware derivative map remains separate.

All-higher-loop and all-energy errors, finite-gravity IR/Regge
and Delta, absolute coupled renormalization, corrected
cosmological cones and a common controlled bounce parent
remain open. No finite-order result is a UV completion.

Scoped P8(a) and A.20-A.23 are unchanged.
Original P8(b), and therefore original P8, is not finished or closed.
