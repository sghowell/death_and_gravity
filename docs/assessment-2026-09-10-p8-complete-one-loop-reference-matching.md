# P8: complete scalar-fermion one-loop canonical matching

Original P8 remains OPEN. This follows the
[first fermion four-scalar audit](assessment-2026-09-10-p8-fermion-four-scalar-increment.md).
No user intervention is needed for the current research.

## One explicit reference scheme for both sectors

[S6.133](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/scheme/FORMULATION.md)
combines the complete old scalar one-loop amplitude and the first
fermion four-scalar contribution at the named MS interaction boundary.
The light mass, residue, H one-point and vacuum energy have explicit
physical local-reference conditions. This is not pure MS for every
relevant parameter, and the internal free scalar still has mass one.

The old scalar subtraction uses I0=Ibar_MS+ell,
ell=log(mF^2)=400 log(10), and the previously fixed finite
quartic contact sigma. Its total counterterms are distinct from
bare coupling counterterms once the field normalization is included.

Let g=G^2, r=Pi_scalar'(1)>0 and fp=f_fermion'(1)>0.
Independent bare-field and bare-coupling equalities give the
canonical-star interaction parameters through first loop order:

    delta G_star=-L G ell/(2Q)+(fp-r)G,
    delta M_star=-g ell/(2Q),
    delta L_star=-3L^2 ell/(2Q)-sigma+2(fp-r)L.

A separate variation of the total counterterms reproduces the
same amplitude conversion. With C(z)=-L+g/(M-z),

    A_new=A0+A_scalar,old^(1)+A_fermion^(1)
          +ell sum_z C(z)^2/(2Q)+sigma+2(fp-r)A0

through one loop. The already canonical vertices receive no
second LSZ factor. The same ultraviolet field/coupling ledger
reproduces the inherited fermion contribution to the MS quartic flow.

The scalar and fermion one-loop graph ownership is complete here:
there is no direct H Yukawa or Phi gauge vertex supplying an
additional one-loop graph. Changes inside these one-loop graphs
caused by the finite parameter conversion are later-order work,
not silently included in this result.

## All local one-loop references are retained

The H one-point counterterm cancels the stationary-heavy tadpole
piece. The remaining scalar two-point increment is

    f_scalar(s)=L T_light/2-g B_HL(s)=-Pi_scalar(s).

Its sign is opposite to the named fermion inverse increment.
The finite light reference mass is 1+Pi_scalar(1)-f_fermion(1);
its large negative numerical value is not a physical tachyon
or a negative free propagator mass.

The full kinetic normalization is kappa=1+r-fp. After the
complete local on-shell subtractions,

    Gamma_can(s)=1-s+[f_fermion,R(s)-Pi_scalar,R(s)]/kappa.

The report also retains both real-scalar vacuum constants and
all fourteen Dirac flavors' constants in the one-loop energy
reference. Massless gauge/ghost one-loop vacuum integrals are
scaleless in the stated regulator; no claim about nonperturbative
gauge vacuum energy follows.

## Quantitative complete one-loop result

A positive-series rational enclosure fixes the actual logarithm,
900<ell_lower<ell_upper<1000. No order-one scale guess is used.
For D=M-2 the scale-conversion coefficient relative to the tree
b2 is -2L+3g/D=g(4-3D)/D^2<0, so the interval arithmetic reverses
the appropriate log and Q endpoints.

The full scalar bound, fermion-box bound and both field slopes
give a complete formal one-loop forward coefficient within
10^-6 relative error of its strictly positive tree reference.
The new conversion and fermion relative contribution is below
10^-202; the inherited scalar bound dominates numerically.

On the radius-two light-pole disc, the combined remainder
coefficient is below 2 times 10^-405. It fixes a mass-one,
unit-residue pole and excludes other poles on the unit disc;
the local Phi potential curvature is positive. These are
one-loop statements, not a full quantum-potential or
all-energy spectrum theorem.

## Independent verification

The native report pins 18 source/proof/test files and 20 fields:
57 named identities and scalar entries, 36 proof gates,
9 controls and 61 rejected inputs.

Final private science: 195 tests in 22.16 seconds.
Independent repository science: 195 tests in 22.32 seconds.
Ordinary replay: 220 tests in 2148.09 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 18892 tests in 2945.99 seconds.

All 519 captured test files are present and unchanged.
Path-list SHA-256:

    99f261487d0691643fc62d54d0474bbe0c0ee6d92c5681fcbc1450fe492046c4

The adapter passed 128 original tuple comparisons. Full-run counters:
33971 domain fallbacks, 7340 exact descents, 4 mixed
fallbacks. The captured snapshot predates S6.134's test files.

Native, ordinary, CLI and both direct science runs use unmodified
SymPy. Only full regression uses the separately audited exact GCD
adapter. The [replay wrapper](p8-replay-runtime.md) supplies
interpreter-only recursion and trusted exact-integer formatting
allowances without changing frozen scientific inputs.

The 101388-character report was transferred losslessly and all
18 source hashes verified. Report SHA-256:

    39613ad0083b068662d7b14475aaf3054db2631dfe142f4ca81343398701805b

Exact 22-file staging contains only this checkpoint, this root
audit, CLAIMS.md and README.md, with staged source/report hashes
verified. Nested continuations and unrelated P4/P9 work are
excluded. The analytic arguments remain written proofs,
not proof-assistant formalization or peer review.

## Active continuation and remaining original obligations

S6.134 extends the first fermion insertion bound to an unbounded
complex half-plane and an explicit finite Euclidean reference
window. S6.135 integrates one precisely named spectral outer-bubble
family after its local subtractions. Both checkpoints are frozen
and native-certified, with independent direct science passes and
fresh ordinary, CLI and full-snapshot replays running.

Neither checkpoint transfers the old pure-scalar two-loop error
budget to the enlarged model. Other two-loop families, finite
reference conversions and later-loop/truncation errors remain.
V also needs its justified high-energy contour and light-cut
treatment; G needs a regulated quantitative Regge/IR allowance;
B needs a common-parent field, derivative, state and cutoff
dictionary with errors.

These are the adopted finite-EFT and necessary-positivity
obligations, not an added demand to construct an all-orders UV
theory. Scoped P8(a) and A.20-A.23 stand unchanged.
Original P8(b) and P8 remain open.
