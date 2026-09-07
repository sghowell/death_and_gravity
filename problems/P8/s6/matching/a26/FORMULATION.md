# S6.19.A26: the literal June-2026 constant-clock domain

## Result and evidence boundary

This checkpoint concerns the **literal coefficient continuation** in
[An et al., arXiv:2606.03302v1](https://arxiv.org/html/2606.03302v1), not a
separately repaired action. It proves a local, integration-by-parts-invariant
failure of a generic constant-clock quadratic linearization. It does not
prove a ghost, a cosmological instability, or a universal UV obstruction.

Use the source signature `(-,+,+,+)` and `X=g^{mu nu} phi_mu phi_nu`, without
a factor of one half. The source bounce has `phi=t, X=-1`. A Lorentz-invariant
constant-clock vacuum would instead have `X=0`.

At a finite field value where the displayed coefficient functions are
defined, put

\[
 f=F_2(\phi_0,0)=\tfrac12-g_1(\phi_0),\qquad
 r=2g_1(\phi_0)-a_1(\phi_0).
\]

For `f != 0`, the exact pole residues are `(r,-r,-r^2/(2f))` for
`(A3,A4,A5)`. On a fixed Minkowski metric, with
`phi=phi0+epsilon psi` and strictly timelike `d psi`, the leading
higher-derivative term is

\[
 \epsilon^2 r\,\mathcal L[\psi],\qquad
 \mathcal L=\frac{L_3[\psi]-L_4[\psi]}{X_\psi}.
\]

For `r != 0`, this functional is not a quadratic form, even modulo
boundary terms. Its Euler expression on the checked two-spatial-direction
quadratic jets is nonlinear in `psi`. A separate compact-time Fourier/IBP
calculation has a nonzero fourth finite difference. These are full-action
derivative-sector tests, not an inference from a pole in one coefficient.
The `A5` pole contributes at order `epsilon^4`, not `epsilon^2`.

At smooth ordinary-source coefficient points, adjoining `F` and the
interacting-chi terms cannot repair this vacuum Hessian: their leading
terms are ordinary quadratic forms with fewer derivatives. The literal
fractional source coefficients at `phi=0` are handled separately below;
no smooth full-potential Taylor remainder is assumed there. The exact
restriction `chi=B(phi)` cancels their reconstructed contribution off shell
and gives the same obstruction at `chi0=B(0)`.

For the named benchmark `epsilon_background=10, w=2, u=1/10`, an analytic
numerator controlling `r` has derivative `-1/20` at zero. Thus `r` cannot
vanish on an open field interval. The theorem therefore excludes a smooth
full-action extension that is unchanged at nonzero `X` near zero on **any
open finite-phi neighborhood**. It does not exclude a pointwise exception
at an isolated residue zero, and makes no nondegenerate claim at `f=0`.

## Literal source choices and cancellation controls

The v1 printed `A4` parenthesis differs from the earlier Ia formula. The
code transcribes both. Their difference is

\[
 \frac{3F_2(1-A_1^2)}{2(F_2-XA_1)^2},
\]

which is regular in `X` when `f != 0` and does not change this theorem.
No cosmological degeneracy calculation is silently repaired or endorsed.

The optional regularity audit explicitly chooses `b=0.018=9/500` in the
printed Q formula. This value is not attributed as a printed numerical
choice: it gives `n_s=24/25` at `epsilon_background=10` via source (21).
Here `Q` is `C2` but not `C3` at `phi=0`. With the named nonconstant positive
entropy background `B`, `W2` is `C1` but not `C2`. These are statements
about literal coefficient regularity in this coordinate presentation.
The entire reconstructed Q-dependent action cancels for `chi=B(phi)`;
after the smooth shift `xi=chi-B(phi)`, its linear source also cancels on
the background. Thus these regularity facts are **not** certified physical
instabilities or a calculation of the reduced perturbation spectrum.

## Matching and verification contract

The frozen adopted S6 formulation and S6.1 smooth-extension certificate
are hash-pinned and replayed read-only. S6.1 permits separately specified
smooth modifications away from the clock tube. Such a modification near
`X=0` changes the literal off-tube action and is not excluded here. We do
not require global real analyticity of a prospective parent.

The source's interacting, phi-dependent chi sector is not the frozen free
canonical spectator M1. Its cubic-action/bispectrum strong-coupling estimates
are not an established common-parent dispersion, finite-gravity, loop or
UV certificate. This checkpoint does not reconstruct quartic scattering,
match the original C/D rows, close S6, or close P8(b).

The report replays symbolic identities, separate exact Fraction Laurent and
Fourier calculations, strict domain/exclusion controls, and independently
authored covariant audits. The written argument, not a proof assistant,
supplies the variational and analytic continuation reasoning.
