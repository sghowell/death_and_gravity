# Source, normalization and proof boundary

The parent action is the frozen S6.20 VARIABLE action, not a new smooth
vacuum extension. All path references below are relative to this child.
The exact Fourier optimization is derived in [optimization.md](optimization.md);
it does not invoke a cosmological positivity or massive-gravity theorem.

## Physical action and actual source

`../../../FORMULATION.md` and `../../../notes/proof.md` fix the positive
equal Einstein terms, canonical interacting clock, free chi, physical g
metric and local background. Signature is +---, with
`R_B=-6(Hdot+2H^2)` and Einstein action `-M^2 R_B/2`.

The frozen S6.21 source module at
`../../src/p8_variable_response/source.py` derives the conserved g-only
probe and its action normalization. Its `sigma=tau^2 Pi/M^2` and actual
canonical light/heavy source projections are retained here. The endpoint
and fundamental matrices are those of the full physical two-tensor
problem; a source coupled directly to only the canonical light field
would be a different calculation.

## Prepared-sector and finite-parameter input

The frozen S6.23 proof at `../notes/proof.md`, Sections 5-6, provides the
actual endpoint inclusion and fixed-source estimates, including
`||H_delta-L||<=200delta`, `||T_delta||<=42`, the `8600delta` regular-light
bound and the `126000000delta ||sigma||1` source term. This child uses
the same punctured interval and K in [1,4]. The new factor 12600000 is
only `sqrt(|J|)=1/10`; it is not a changed source coupling.

The comparison with the exact analytic prepared output has coefficient
8400 because it uses `T_delta H_delta^-=H_delta^+` directly. The Lv and
H_delta^+v targets are explicitly distinguished. The delta=0 equations
used for loading are restricted to J; no physical delta=0 action is
asserted through u=0.

S6.26's actual causal response and S6.27's fixed-positive-pulse phase result
remain separate frozen source-aware statements. A signed optimized
preparation is not placed in S6.27's positive pulse class, and neither
result is overwritten by a projected-source prescription.

## Mathematics established in this child

The actual Gramian and H0^2 right inverse have independent action/ODE
proofs and exact numerical/arithmetic replays in the companion files.
Their constants are inputs to the Fourier arithmetic, not premises
certified merely by entering numbers in `spectral.minimum_cost_envelope`.

Plancherel, bounded-band Fubini, compactness of the continuous interval
kernel, and the elementary entire-Fourier-transform argument yield the
positive tail operator. The optimization note gives the Hilbert-space
projection, multiplier/dual, strict-budget density/correction and endpoint
attainment arguments explicitly. The all-order Hilbert-space theorem is
not claimed to follow from the finite-dimensional KKT regression alone.
Likewise the continuous sine Taylor remainder, not sampled frequencies,
justifies the sinc error bound. No external source supplies a computed
optimal band-limited physical preparation for this bimetric branch.

The validated midpoint ODE uses python-flint 0.9.0 with FLINT 3.6.0,
256-bit Arb balls and an explicitly proved interval Taylor remainder.
The numerical trust boundary is the documented outward ball arithmetic,
not a floating-point integrator's error estimate. The relevant primary
API references are [Arb real balls](https://python-flint.readthedocs.io/en/latest/arb.html),
[Arb matrices](https://python-flint.readthedocs.io/en/latest/arb_mat.html) and
[Arb power series](https://python-flint.readthedocs.io/en/latest/arb_series.html).
The final matrix's positive Schur pivots are also checked independently
using exact Fraction interval arithmetic after rational outward rounding.

Actual loading moments for the band Gramian and a numerical optimizer
are left uncomputed. The certificate is intended to replay their exact
finite-rank formula and continuous error, together with proved cost
envelopes and constructive controls, not to manufacture those missing
numerical results. The distinction between an attained L2 optimizer and
an H0^2 infimum is part of the claim, not a software limitation.

No scattering-positivity theorem is applied. The adopted S6 contract at
`../../../../../FORMULATION.md` still requires separate vacuum, matching,
interaction, finite-gravity and remainder evidence. A physical source
Fourier-energy bound is not a rolling mode gap, cutoff or UV completion.
