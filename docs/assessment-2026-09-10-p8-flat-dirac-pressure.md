# P8 continuation: complete specified homogeneous flat Dirac tensor

S6.169 completes the pressure component of the specified free flat
one-loop tensor. Original P8(b), V/G/B and original P8 remain OPEN.

## Result and reference boundary

The physical pressure operator has both a diagonal energy component and
an off-diagonal mass-coherence component. Four exact rotating frames
control both. Their complete half-line state/projector remainder,
integrated over all momenta, is below1e410.

The full d=3-2epsilon angular measure and isotropic factor1/d are kept
until covariant counterterms cancel the poles. At Q=16pi^2,
ell=log(M^2/m^2), the finite two-derivative terms are

    rho2=-Mdot^2(ell+2/3)/Q,
    P2=[2M Mddot ell/3-Mdot^2 ell/3+2Mdot^2/3]/Q,
    P2-rho2=partial_t^2{M^2[log(M^2/m^2)-1]}/(3Q).

The last term is a curvature improvement: evaluating R=0 before metric
variation would incorrectly discard it. Additional finite curvature
terms are explicitly set to zero in this common flat MS prescription;
flat stress does not determine a complete curved gravitational reference.

The derivative pressure is below1e595. With the same entire potential
and fixed saturated mass reference from S6.168, at every real time

    |P|<1e789,  |rho|<1e789,
    max(|rho|,|P|)/kappa<1e-11.

These are all components of the homogeneous tensor in its fixed
orthonormal rest frame, not a bound uniform under arbitrary boosts.
The ratio to kappa is not a relative error against zero bounce density.
This is a quadratic/reference contribution in the MS scalar coordinate,
not the full interacting/canonically re-expressed parent tensor.

## Verification

Final private science passed193 tests in22.16 seconds.
Fresh repository science passed193 in21.06 seconds.
Ordinary replay passed218 in2024.42 seconds; independent native CLI
replay passed. The complete captured P8 snapshot passed28470 tests
in3281.82 seconds, final exit code0.

All591 captured test files remained present and unchanged.
Path-list SHA-256:

    f0a2de5c5ae3fa990cf6d7ee1f6907e25d2fb42ae7815435243e90b2889fed51

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34093 domain fallbacks,7340 exact descents,88 mixed
fallbacks. Native, direct science, ordinary and CLI use original SymPy.

The native29415-character report was transferred losslessly in three
chunks and independently checked against15 source/proof/test hashes.
Its20 fields record37 named identities,37 scalar entries,22 proof gates,
9 controls and93 rejected inputs. Report SHA-256:

    97f1821d03f6e2f60367caf25bbc650b7122e31aa07cb77409755c084b89a972

Independent tests cover full-function Pauli projectors, finite angles,
600-digit coherence cancellation, complete radial integrals,
nonzero-regulator angular/Gamma factors, full mass-profile derivatives
and finite-improvement negative controls. They supplement written
analytic estimates, not validated numerical integration.

An initial private scope-wording assertion was corrected before freeze.
No numerical bound was weakened and no prior frozen scientific byte
was changed. Exact19-file staging excludes later continuations and the
unrelated P4/P9 changes.

## Next frontier

The actual CD curved in/out Hadamard-state checkpoint S6.170 has passed
native and direct science checks; its fresh replay set is running.
A separately specified absolute curved one-loop tensor, including the
evanescent curvature prescription and Newton reference, is in private
development. Neither result by itself supplies the complete interacting
parent state, cutoff, quantum target matching or background/response.

Vacuum contour/truncation, finite-gravity IR/Regge estimates and the
remaining common-parent matching obligations are research work, not a
user-choice blocker. Scoped P8(a) and A.20-A.23 are unchanged.
