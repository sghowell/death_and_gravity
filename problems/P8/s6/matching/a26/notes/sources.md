# Primary source and convention audit

The local variational and Fourier proofs in `notes/proof.md` are new
calculations in this repository. The papers below supply the literal
input action, not the conclusion that its constant-clock chart fails.

## Literal June model

[An, Hong, Kang and Mun, arXiv:2606.03302v1, 2 June 2026](https://arxiv.org/html/2606.03302v1):

- Section 1 fixes `(-,+,+,+)` and reduced Planck units.
- Section 2.1, equations (1)–(4), supplies the complete two-field action,
  elementary Li invariants, and exceptional relations. `X` has no half.
- Section 2.2, (17), (19), (23)–(24), (28), (30)–(31), supplies the clock
  background, reconstructed source, coefficient profiles and benchmark.
  `cosh^{-1}` in B denotes its reciprocal, not inverse hyperbolic cosine.
- Section 3.1–3.3 studies cosmological perturbations and cubic/bispectrum
  strong-coupling estimates. These are not a flat-vacuum positivity theorem
  or a controlled gravitational decoupling/loop remainder.

The source has a phi-dependent interacting chi sector, not the frozen
canonical free M1. Its background `X=-1` is not the constant-field `X=0`
configuration tested here. The fractional tilt audit explicitly chooses
`b=0.018`, which source (21) maps to `n_s=0.96` at epsilon=10; neither
number is claimed as a printed numerical choice in the paper. The
nonconstant B example is from (35)–(36). No divergence of W2 is inferred.

## Printed A4 discrepancy, explicitly retained

[An et al., arXiv:2501.09985v2](https://arxiv.org/html/2501.09985v2),
section 2.1 equation (3b), has the term
`4(3F2+16 X F2X) A1^2`. The June v1 HTML (4b), also checked against its
[PDF](https://arxiv.org/pdf/2606.03302), instead prints
`4(3F2+16 X F2X A1^2)`. The code's `printed=True` is the latter.
Their symbolic difference is regular near the nondegenerate X=0 point
and does not affect any leading-domain conclusion. We use the earlier
formula only as an explicit comparison; we do not silently repair the
June action or transplant the earlier paper's opposite metric signature.

## Frozen repository context

Paths here are relative to the A26 child root.

- `../../FORMULATION.md`: adopted S6 common-parent/vacuum/physical-matter
  contract, SHA-256
  `d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901`.
- `../../certificates/vacuum-extension.json`: S6.1 smooth off-tube extension
  context, SHA-256
  `6a7ab2a7bf6719d3829795779a876dab0694a5b99383726c14f0880515d2d7d4`.
  Its own read-only verifier and prior arithmetic lineage are replayed.
- `../../notes/vacuum-extension.md`: the distinction between rolling-tube
  data and a separately specified smooth vacuum extension. It is already
  protected by the S6.1 source manifest, not edited or newly assumed here.

No dispersion theorem, positivity sign, observational inference, or UV
noncompletion claim is imported from secondary summaries.
