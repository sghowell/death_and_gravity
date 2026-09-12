# P8 continuation: full constrained ADM Proca vertices

S6.214 establishes local and finite-regulator inputs for lapse, shift and arbitrary spatial metric directions. It does not construct their full continuum current or close original P8.

**Current correction:** the [retarded phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) withdraws the physical S212/S213 finite-matching/current identifications. The new S214 ADM algebra below is independently derived from the Hamiltonian and is unaffected by that endpoint-branch error. Its frozen inherited status table is a historical record, not renewed certification of those identifications.

## Result and scope

For the original +--- metric, use lapse N, shift beta and h=a^2 exp(Q). Retaining the temporal Proca constraint gives

    A0 = beta.A - N div(Pi)/(m^2 sqrt(h)),

and the Hamiltonian includes the full positive longitudinal constraint term and

    H_shift = beta^i [Pi^j F_ij - A_i div(Pi)].

The latter is equivalent, modulo the stated spatial divergence, to Pi.Lie_beta A. Dropping its longitudinal transport term changes the physical shift vertex.

At three spatial dimensions write B_Q=Q-tr(Q)I/2. The exact lapse Hamiltonian is

    H_N=N/2 [Pi.exp(B_Q).Pi/a
             +curl(A).exp(B_Q).curl(A)/a
             +a m^2 A.exp(-B_Q).A
             +exp(-tr(Q)/2) div(Pi)^2/(a^3 m^2)].

Every first and mixed second ADM vertex is retained, including noncommuting ordered exponential words, trace-constraint terms and lapse/spatial cross contacts. The lapse is linear; the independent ADM chart has no lapse/lapse or shift-containing second Hamiltonian vertex.

That last fact does not remove the nonlinear four-metric chart contact. The physical metric map has

    h_D00=2n_D, h_D0i=-a^2 beta_Di, h_Dij=-a^2 Q_Dij,
    q_DG00=2n_D n_G-2a^2 beta_D.beta_G,
    q_DG0i=-a^2(Q_D beta_G+Q_G beta_D)_i,
    q_DGij=-a^2(Q_D Q_G+Q_G Q_D)_ij/2.

The original nonzero one-point current contracts with all appropriate chart terms.

The complex Fourier shift form has exact norm a|beta|. With

    sigma(D)=|n_D|+(5/2)||Q_D||op+a|beta_D|,

the relative first and second Hamiltonian forms are bounded by sigma(D) and sigma(D)sigma(G), uniformly in both internal momenta. Complex Fourier blocks are not presumed individually Hermitian. Conservative finite-band memory/contact bounds still grow as K^5 and K^4 and are not ultraviolet limits, renormalized Ward identities, reduced canonical norms or inverse estimates.

## Verification

Initial core ADM, vertex, chart and norm checks passed. Independent science first passed41 checks, then was expanded to62. One expanded test compared an unexpanded complex reversed-pair expression structurally and failed despite all36 exact residuals vanishing. The test was corrected before freezing to use exact cancellation; no physical formula was changed.

Integrated science passed240 tests in6.80 seconds; preflight240 in2.62 seconds; final private240 in6.64 seconds; repository240 in6.89 seconds. Lint passed.

Fresh original-SymPy ordinary replay passed265 tests in2299.56 seconds; CLI passed. The complete captured P8 regression passed42063 tests in4331.58 seconds with final exit code0. It retained681 captured files, path-list SHA

    d91acfdfa0ce8e3ca578c0c27a37329b8e07de015d14995879cfa9938966f8c1

The full-only exact-GCD adapter passed128 comparisons with original tuple results. Final counters were35531 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct science, ordinary and CLI used original SymPy.

The37581-character native report was transferred in four verified chunks. Its twenty fields contain39 named identities,697 scalar entries,35 gates,nine controls and136 rejected inputs. All eighteen frozen source hashes and the report SHA were checked:

    2d2600e003055aec80398fdafa69e8e01a54c7c52e915832a51a0d23aa12db01

Publication contains exactly22 files and excludes S215, private S216 and unrelated P4/P9 work. Continuous analytic estimates are written proofs, not FORMALIZED. The earlier tests' passing status is not offered as a cure for the phase error; the independent local derivation is the relevant evidence for this claim.

## Continuing frontier

The conditional prepared Ward reconstruction in frozen S215 has passed ordinary and CLI checks; its complete regression remains pending at this audit. Its previously claimed known tracefree input requires the explicit S213 correction.

Private S216 has independently diagnosed the branch error, reversed odd endpoints consistently, matched all three ordered scalar physical poles, and verified48 full-volume curvature Hessians at spatial dimensions3,4,5,6. Full scalar finite matching, corrected current assembly, scalar anchor and state/time/contact bounds remain under review.

Original P8(a)'s scoped closure and A.20-A.23 are unchanged. Full reduced response/inverse, controlled nonlinear quantum background, remaining parent/loop/heavy/cutoff matching, vacuum cuts/contour/truncation and finite-gravity IR/Regge remain open. No user-intervention blocker has been identified.
