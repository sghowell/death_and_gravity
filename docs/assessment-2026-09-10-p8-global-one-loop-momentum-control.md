# P8: global internal-momentum control for the one-loop kernel

Original P8 remains OPEN. This follows the
[complete one-loop matching audit](assessment-2026-09-10-p8-complete-one-loop-reference-matching.md).
No user intervention is needed for the current research.

## A bound that covers internal loop momenta

[S6.134](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/scheme/internal_momentum/FORMULATION.md)
starts from the same complete fermion two-point kernel and physical
on-shell references. Its new bound is not restricted to small external
momentum compared with the fermion mass.

With A=x(1-x), the parameter kernel has the exact decomposition

    g_x''(s)=-A/(mF^2-A s)
             -mF^2 A(1-4A)/(mF^2-A s)^2.

Both numerator weights are nonnegative; their integrals are 1/6
and 1/30. For Re(s)<=sigma, 1<=sigma<4mF^2, the parameter
gap d_sigma=mF^2-sigma/4 is positive. Subtract the whole
affine regulator reference and both finite on-shell anchors first.
Taylor's exact integral remainder gives

    |f_R(s)/(s-1)^2|
      <= (NY/Q)[1/(6d_sigma)+mF^2/(30d_sigma^2)].

The quotient has a removable extension at one. The estimate is
uniform on an unbounded complex half-plane, without a momentum
expansion. At sigma=1 its canonically normalized actual upper
bound is approximately 5.36883318866 times 10^-608.

The selected symmetric outer-bubble routing satisfies

    Re(s_line)=-(q0 minus/plus Im E)^2-|q_spatial|^2+(Re E)^2
              <=3/4<1

for all real q on the unit forward disc, E=sqrt(s)/2.
This is a pointwise identity, not a complex contour translation.
Other graph routings are not automatically covered.

## Retaining large-momentum decay and both inverse signs

For s=-t, d=t+1 and K=4mF^2-1, the same exact subtraction yields

    0<=-f_R(-t)<=(2NY/Q)d log(1+d/K).

Dividing by d^2 retains a logarithmic-over-momentum decay envelope
for the inserted free propagator. A coarser constant bound is
also valid but cannot by itself justify a complete outer integral.

The scalar remainder obeys 0<=Pi_scalar,R(-t)<=d r.
Combine both sectors with kappa=1+r-fp, keeping their distinct
inverse signs. The selected complete one-loop inverse has
fractional defect at most

    [r+(2NY/Q)log(1+(t+1)/K)]/kappa.

On the explicitly named test range 0<=t<=(10^400)^2, with
mF=10^200, the logarithm is bounded by the inherited exact
400 log(10) enclosure. The actual rational defect upper bound
is approximately 4.94514301641 times 10^-204, below 10^-202.
Thus this selected one-loop Euclidean inverse has no zero
throughout that finite range.

This mathematical reference-energy window is not a derived
Wilsonian cutoff, an all-order result, reflection positivity,
or the full Lorentzian/gauge spectrum.

## Independent verification

The native report pins 18 source/proof/test files and 20 fields:
51 named identities and scalar entries, 31 proof gates,
9 controls and 196 rejected inputs.

Final private science: 384 tests in 22.02 seconds.
Independent repository science: 384 tests in 22.75 seconds.
Ordinary replay: 409 tests in 2210.39 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 19301 tests in 2915.85 seconds.

All 521 captured test files are present and unchanged.
Path-list SHA-256:

    1d7865ab47ca353c1906e9b1878ab112bc44c92293b9b0d4980e5130b73b79ab

The adapter passed 128 original tuple comparisons. Full-run counters:
33992 domain fallbacks, 7340 exact descents, 4 mixed
fallbacks. The snapshot predates S6.135's test files.

Only full regression uses the audited exact GCD adapter.
All native, ordinary, CLI and direct scientific runs use
unmodified SymPy with the documented interpreter allowances.
An initially truncated export was discarded and re-exported;
no incomplete report was saved.

The final 48336-character report was transferred losslessly and
all 18 source hashes verified. Report SHA-256:

    a6d91c84be131cfd36d877a870213f15c61c16e487ec00816c1eff98cfb98404

Exact 22-file staging verifies staged source/report bytes and
contains only this checkpoint, this root audit, CLAIMS.md and
README.md. It excludes nested continuations and unrelated P4/P9
work. The analytic arguments are written proofs, not proof-
assistant formalization or independent peer review.

## Active continuation and remaining work

S6.135 derives a positive spectral representation of the same
insertion and bounds a local-quartic outer-bubble contribution.
S6.136 extends it to the once-differentiated scalar quadratic
sector, including local terms and the mixed heavy-light bubble.
Both are frozen and native-certified, with passing independent
science and fresh replays underway.

A new private two-loop ownership analysis organizes the
remaining fermionic primitive terms and their separate local
counterterm/finite-reference conversions. A complete inventory
is not a complete error bound, and the existing subset results
are not being promoted to the full two-loop amplitude.

The adopted V contour and error obligations, finite-gravity
IR/Regge allowance for G, and common bounce-parent field/state/
cutoff control for B remain open. Scoped P8(a) and A.20-A.23
are unchanged. This is still the agreed finite-EFT/necessary-
positivity task, not an added requirement to construct an
all-orders UV theory. Original P8(b) and P8 remain OPEN.
