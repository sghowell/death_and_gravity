# P8: first fermion spectral outer-bubble contribution

Original P8 remains OPEN. This follows the
[global one-loop momentum audit](assessment-2026-09-10-p8-global-one-loop-momentum-control.md).
No user intervention is needed for the current research.

## From the complete insertion to an outer integral

[S6.135](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/threshold/four_scalar/scheme/internal_momentum/spectral_bubble/FORMULATION.md)
derives a positive spectral representation of the whole first
on-shell-subtracted fermion propagator insertion. With
T=4mF^2, C=2NY/Q, Q=16 pi^2 and beta=sqrt(1-T/u),

    -f_R(s)/(s-1)^2 = integral_T^infinity w(u)/(u-s) du,
    w(u)=C u beta^3/(u-1)^2 >= 0.

The proof starts with the exact parameter kernel, differentiates
the proposed representation twice, and matches both on-shell
anchors. Its analytic continuation and compact majorants are
explicit. No high-energy contour assumption is used.

This is the first formal propagator insertion, not an exact
resummed propagator. Its spectral weight has infinite total
measure: for u>=2T, w(u)>=C/(4u). The weighted propagator
integral converges, but this does not prove a normalized exact
Kallen-Lehmann measure or reflection positivity.

## A renormalized two-loop family and its error

Insert this spectral propagator into the scalar bubble with
two literal local quartic vertices. There are three channels
and two possible internal insertion positions. The two positions
cancel the original bubble symmetry factor one half.

The proper inner on-shell reference is paired first. Remove the
remaining overall quartic constant at channel invariant two
before integrating over the spectral mass. For Re(z)<=4,
u>=T>=16, the resulting parameter denominator is strictly
positive in real part. Its derivatives obey

    |B'(z,u)| <= 1/[2(u-4)],
    |B''(z,u)| <= 1/[3(u-4)^2],
    w(u) <= 4C/u.

These estimates prove convergence of the once-subtracted outer
integral and justify the differentiated expression. The two
crossed forward channels have Taylor factors one half; their
sum leaves one second derivative. The t-channel constant has
zero quadratic coefficient. Thus

    b2_family = (L^2/Q) integral w(u) integral A^2/Delta_2^2 > 0,
    b2_family <= NY L^2/(3 Q^2 mF^4).

The actual rational upper bound is approximately
3.25528355986 times 10^-1419. Relative to the tree coefficient
4 lambda it is below 10^-819. The same estimate controls
the quadratic remainder on the unit forward disc.

Changing the finite overall local quartic contact cannot change
this family's b2. That observation is not extended to nonlocal
heavy-vertex families or their separate finite references.
The calculation has no extra resummed field factor; changes to
its leading couplings start at the next formal loop order.

## Independent verification

The native report pins 18 source/proof/test files and 20 fields:
52 named identities and scalar entries, 27 proof gates,
9 controls and 54 rejected inputs.

Final private science: 244 tests in 21.37 seconds.
Independent repository science: 244 tests in 21.71 seconds.
Ordinary replay: 269 tests in 2330.04 seconds.
Independent native command-line replay: passed.
Complete P8 snapshot: 19570 tests in 2943.42 seconds.

All 523 captured test files are present and unchanged.
Path-list SHA-256:

    30bca9f2eb4dcb0a6e84ec53bac25792530602d03683c3a673b186000d0c175b

The adapter passed 128 original tuple comparisons. Full-run
counters: 33982 domain fallbacks, 7340 exact descents and
4 mixed fallbacks. This snapshot predates S6.136's tests.

Only full regression uses the audited exact GCD adapter.
All native, ordinary, CLI and direct science runs use
unmodified SymPy with interpreter-only runtime allowances.
The 27596-character report was transferred losslessly,
and all 18 source hashes verified. Report SHA-256:

    be4188f2bcd3ff3b0d68a20fd2b385c0c9e936c2ab4df52626cefabf8e21e5c2

Exact 22-file staging checks staged report/source bytes and
excludes nested continuations and unrelated P4/P9 changes.
The analytic proofs are written arguments, not proof-assistant
formalization or independent peer review.

## Continuation and scope

S6.136 has independently passing science and native certification;
its fresh replays continue. It treats the paired quadratic
sector, including the mixed heavy-light bubble and local
reference cancellations. S6.137 now gives a native-certified
complete primitive fermion-sector ownership ledger through
four scalar derivatives, with its fresh replays underway.

These identify and bound selected contributions, not the
complete enlarged-model two-loop amplitude. Remaining primitive
integrals, proper counterterms, finite canonical conversion and
later-loop/truncation control are still research obligations.
The adopted V contour and error conditions, finite-gravity
IR/Regge allowance for G and common bounce-parent field/state/
cutoff dictionary for B remain open. Scoped P8(a) and A.20-A.23
are unchanged. No all-orders UV construction has been added to
the agreed finite-EFT/necessary-positivity contract. Original
P8(b) and P8 remain OPEN.
