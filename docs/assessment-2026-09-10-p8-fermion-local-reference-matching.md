# P8: one-loop fermion local-reference matching

Original P8 remains OPEN. This follows the
[first light-cut subtraction audit](assessment-2026-09-10-p8-first-gauge-light-cut-subtraction.md).
No user intervention is needed for the current research.

## The local thresholds cannot be inferred from the tiny cut

[S6.131](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/FORMULATION.md)
computes the complete one-loop fermion increment of the
prospective GY14-unbroken reference boundary. It retains
both opposite Yukawa masses, all twelve inert Dirac
flavors, three colors, and the whole dimensional
reference before subtraction.

For N=6, Y=y^2, Q=16pi^2 and reference scale nu=mF,

    V_F(0)=63mF^4/Q,
    f(0)=V_F''(0)=4NYmF^2/Q,
    v4=V_F''''(0)=-64NY^2/Q.

The inert flavors do not affect Phi derivatives, but
do affect the fixed vacuum-energy reference. At the
named boundary the mass threshold is between 10^193
and 10^195 and the vacuum constant between 10^799
and 10^800. Neither is bounded by the tiny gauge
cut coefficient. No naturalness exclusion is imposed.

Every even field degree beyond four has positive
coefficient 96Nc/[k(k-1)(k-2)(k-3)(k-4)] on the
fermion-gap domain |y Phi_reference/mF|<1. The
complete higher-field sum has an explicit bound on
the closed half-domain. This is a one-loop field
series, not an all-loop or bounce-field bound.

## On-shell reference and field map

An independent Dirac trace gives the complete two-point
increment, including its tadpole and affine UV reference.
Subtract f(1)+(s-1)f'(1) only after retaining those terms.
Grouped exact moment series enclose the anchors and
the tiny remaining curvature term without subtracting
huge decimal approximations.

Fix the vacuum constant, set m_ref^2=1-f(1), and let
kappa=1-f'(1). The constant map
Phi_can=sqrt(kappa) Phi_reference gives

    Gamma_can(s)=1-s+f_R(s)/kappa.

On |s-1|<=2 the full parameter-integral remainder
obeys |f_R(s)|<=B|s-1|^2. At the actual boundary,
0<f'(1)<10^-205 and B/kappa<10^-606. The selected
kernel therefore has a mass-one unit-residue pole
and no other pole on the unit disc.

The same map changes G,y and the potential quartic.
The fermion quartic shift is below 10^-8 times the
tree heavy-eliminated quartic margin. The heavy
square retains the frozen parent's cubic sign, which
is checked directly. Together with the positive
curvature and higher even terms, this proves positivity
away from the origin of the selected scalar-tree plus
one-loop-fermion potential on its finite gap domain.

This is not the full new-model quantum potential.
Canonical couplings are not unchanged MSbar ray
numbers, the negative local reference mass is not
a physical tachyon, and the large counterterm is not
inserted as a negative free scalar propagator mass.
Old scalar loops and later mixed/gauge loops still
need their own matching and error accounting.

## Independent verification

The native report pins 18 source/proof/test files and
20 fields: 67 named exact identities and scalar entries,
32 proof gates, 9 controls and 30 rejected inputs.

Final private science: 202 tests in 20.32 seconds.
Independent repository science: 202 tests in 20.45 seconds.
Successful ordinary replay: 227 tests in 2150.56 seconds.
Independent native CLI: passed.
Complete P8 snapshot: 18400 tests in 2856.82 seconds.

The first ordinary attempt failed at the default Python
recursion limit during pytest's ancestor import, with
an incomplete collection; it is not counted as a pass.
The read-only retry raised the recursion allowance to
4000 and collected all expected tests. No frozen source
or scientific algorithm was changed. The separately
tested [replay bootstrap](p8-replay-runtime.md) makes
the runtime settings reproducible for later checkpoints.

All 515 full-snapshot test files are present and unchanged.
Path-list SHA-256:

    c426172da9ed73115f11af75e7e06851deead7073bc12473a32296d8d68535c1

The adapter passes 128 original tuple comparisons.
Full-run counters: 33995 domain fallbacks, 7340 exact
descents, 4 mixed fallbacks. The snapshot predates
S6.132's two test files. Only full regression uses
the audited exact GCD adapter; private, repository,
native, ordinary and CLI science use unmodified SymPy.

The 101246-character native report was transferred
losslessly and all 18 source hashes were verified.
Report SHA-256:

    2c8ed16ab546b49992e59dc512050e004b28069c3fee9705a831224bfce0d3a7

Exact 22-file staging includes only the checkpoint,
this root audit, CLAIMS.md and README.md. It verifies
staged native source/report bytes and excludes nested
children and unrelated P4/P9 work. Determinant,
continuum and complex-domain arguments are written
proofs, not proof-assistant formalized or peer reviewed.

## Active continuation and original obligations

S6.132 is frozen and native-certified, with 247 passing
tests in each independent final science run. Its
four-scalar fermion increment retains the local threshold,
complete quadratic-momentum term and every higher
momentum degree, plus the pole-field dictionary.

S6.133 is frozen and native-certified following 195
passing private science tests. It supplies the finite
scheme conversion needed to combine the old scalar
and new fermion one-loop amplitudes. Both checkpoints'
fresh ordinary, CLI and full-snapshot replays are running;
this audit does not promote them early.

The next estimates address internal loop momenta,
where a small external-disc bound is insufficient.
New-model later-loop errors, justified V contours,
finite-gravity IR/Regge allowance and common-parent
bounce field/state/cutoff control remain unproved.
These are the adopted finite-EFT and necessary-
positivity obligations, not a new demand to construct
an all-orders UV theory. Scoped P8(a) and A.20-A.23
stand unchanged. Original P8(b) and P8 remain open.
