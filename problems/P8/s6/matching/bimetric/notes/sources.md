# Primary-source and convention audit

The local calculations are derived from the explicit action in
[FORMULATION](../FORMULATION.md), not imported sign-by-sign from a paper.
The following primary sources provide independent equation and procedure
checks. The certificate does not adopt their broader health or completion
claims as a theorem about CD.

## Cosmological equations

[von Strauss et al., *Cosmological Solutions in Bimetric Gravity and their
Observational Tests*](https://arxiv.org/pdf/1111.1655), equations (2.10),
(2.12), and (2.14)--(2.18), retain the relative metric lapse and couple
matter only to `g`. In their FLRW notation `Y=b`, `X=N_f` after `N_g=1`,
`M_star^2=F/G`, and `m^2 beta1=nu/G`. Their signature is `-+++`, with
mixed matter components `T^0_0=-rho`, `T^i_i=p`. We compare physical
Friedmann equations after this dictionary, not the bare overall action
or stress signs. Their equation (2.12), before the division used in (2.13),
is the appropriate Bianchi relation at a zero of `adot`. Here both lapses
are varied before gauge fixing and the undivided identity is checked
directly, so neither the divided Bianchi chart nor a differentiated
Friedmann equation loses the independent bounce condition.

## Mass eigenstates and correct elimination

[Gording and Schmidt-May, *Ghost-free infinite derivative gravity*](https://arxiv.org/pdf/1807.05011),
equations (2.2)--(2.7), supply a proportional-vacuum spectrum cross-check.
Their normalization maps to `G=2 m_g^2`, `F=2 alpha^2 m_g^2`,
`mathfrak m^4=2 m^2 m_g^2`, with the Einstein sign changed for the P8
curvature/signature convention. At proportionality one, their mass gives
`mu^2=nu(1/G+1/F)`; our determinant and full FP expansion independently
derive it. Section 3.1 distinguishes solving the eliminated field's own
equation from substituting a solution of the other field's equation.
Section 2.3 also explains that `f` is not a pure massive eigenfield, and
warns about strong self-interactions in a heavy small-alpha limit. We
therefore compute true mass-basis integration with its source terms
separately, and do not interpret the Schur spectral disk as a physical
Wilsonian or cosmological cutoff.

## Source terms and equation-versus-action elimination

[Hassan, Schmidt-May and von Strauss, *Higher Derivative Gravity and
Conformal Gravity From Bimetric and Partially Massless Bimetric Theory*](https://arxiv.org/pdf/1303.6940),
section 2.1, equations (2.2)--(2.4), exposes the functional chain-rule term
when a metric is eliminated. Appendix A, in particular (A.7)--(A.10),
retains sources in a coupled scalar example and demonstrates why a
source-source term cannot be discarded when comparing amplitudes. Our
two-metric vacuum calculation is independent of that scalar example:
it uses the actual FP mass, all ten relative metric equations, conserved
spin projectors, and a fixed physical source. The wrong-field substitution
and omitted source functional are tested negative controls, not matching
prescriptions. The primary paper's higher-curvature construction is not
assumed to identify an interacting heavy-mode EFT for CD.

## What is new here and what remains open

The bounce contradiction, its finite-window robustness, all local
normalization factors, source projectors, curvature contractions, and the
stated rational remainder are the explicit derivations certified in this
subtree. None of the sources supplies a regular CD branch for this exact
family, a bound on its rolling heavy modes or quantum remainder, or a
general UV verdict. The zero flat quadratic tree `c_R` is not a claim
about loop counterterms or the full nonlinear covariant effective action.
