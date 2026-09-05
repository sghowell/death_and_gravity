# Why this conditional UV framework

Research checked 2026-09-04. Primary papers below guide applicability, not a
new positivity prior on the old row classification. Our recommendation is
the synthesis in `../FORMULATION.md`; no paper proves that synthesis for P8.

## Vacuum versus background

- [Xu and Zhou, 2306.06639v2](https://arxiv.org/html/2306.06639v2), introduction,
  sections 2, 4 and 5. Their fully crossing bounds address scalar scattering
  with flat-space/decoupling hypotheses. They explicitly caution against
  transferring them to time-dependent observables. The combined Horndeski/
  DHOST discussion also matters: exclusions of restricted DHOST sectors must
  not be transplanted to our arbitrary baseline F. We use neither their
  coefficient dictionary nor their numerical bounds without a new derivation.
- [de Boe et al., 2403.13096v2](https://arxiv.org/html/2403.13096v2), introduction
  and section II. Their phenomenological application explicitly assumes a
  promotion of Minkowski bounds to a cosmological setup. That is a useful
  conditional comparison, but not the certificate-level bridge P8 requires.

## Broken boosts are not the only issue

- [Grall and Melville, 2102.05683v1](https://arxiv.org/pdf/2102.05683v1),
  sections I--II. The boost-breaking construction retains spatial and time
  translations and states crossing, unitarity, analyticity and boundedness
  hypotheses. A rolling FLRW background is not automatically that problem.
- [Creminelli et al., 2312.08441v2](https://arxiv.org/html/2312.08441v2),
  sections 2--3 and conclusions. A weakly coupled complex-scalar model in a
  finite-charge state exhibits momentum-dependent mode mixing and
  nonstandard amplitude analyticity despite its healthy microscopic origin.
  This is a concrete reason not to infer the boost-breaking dispersion
  assumptions from causal low-energy propagation alone. It is not a no-go
  theorem against every boost-breaking positivity construction.
- [Lee and Melville, 2512.20706v2](https://arxiv.org/html/2512.20706v2),
  sections I--III. This newer work constructs positivity tests for in-in
  propagators on fixed de Sitter, retaining background-induced dissipative
  effects. It supplies a useful alternative direction but not a ready-made
  four-point bound for our nonstationary bouncing scalar-tensor geometry.
  Its background, observable and state must not be replaced silently.

## Gravity and infrared sensitivity

- [Tokuda, Aoki and Hirano, 2007.15009v2](https://arxiv.org/pdf/2007.15009v2),
  sections 3.1--3.4. The fixed-negative-transfer analysis combines the
  graviton pole with high-energy Regge contributions. Its finite residual can
  invalidate strict forward positivity. Their approximate size estimate
  additionally assumes scales for Regge residue/trajectory derivatives.
  Hence our contract keeps these assumptions and the finite correction
  explicit, rather than adopting a bare positive pole-subtracted coefficient.
- [Caron-Huot and Tokuda, 2406.07606v2](https://arxiv.org/html/2406.07606v2),
  introduction, section 2.1 and conclusions. Their string-inspired loop
  calculation shows that forward high-energy scattering retains sensitivity
  to light masses; Regge data can cancel infrared-enhanced negative low-energy
  terms. This prevents us from declaring the gravitational error universally
  order one in units of the Regge scale. Their example does not supply a
  numerical correction for P8. We choose the conditional Regge/dispersion
  route, with a stated remainder, as the first finite-gravity test; smeared
  sum rules are a possible later cross-check, not a result implemented here.

## Recommendation and rejected shortcuts

Freeze a vacuum-based, Lorentz-invariant UV target with explicit EFT matching,
then retain gravitational corrections when returning to finite M. This is a
conditional admissibility class, not a claim to cover every possible UV
theory. It is informative only when the vacuum is genuinely related to the
bounce. The new splice example tests this point constructively.

Neither a chosen healthy off-tube vacuum nor a very large M*tau proves the
required connection. Conversely, the failure of the original polynomial
extension to possess a constant-clock vacuum is not a no-go against the
clock tube. Assuming exact global coefficient analyticity would change that
question; its precise conditional obstruction is recorded separately.

No unqualified pointwise bounce inequality, hard-channel-to-forward limit,
P9 sign convention, or fully UV-complete bounce is inferred from these papers.
