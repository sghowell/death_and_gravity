# Source and normalization audit

The physical model and full operator are the frozen sources in
`../../FORMULATION.md`, `../../src/p8_variable_beta/canonical.py`, and
`../src/p8_variable_response/operator.py`, with paths stated relative to
this child root. The directly pinned S6.21 report is
`../certificates/matched-physical-response.json`, SHA256
`7ac41e1ca23c9fe9649e6fc703c5d3d236fb599481f921db319416e10865f118`.
Its read-only verifier replays VARIABLE and the adopted S6 contract.
No frozen source, sign, action, physical metric or source is modified.

The motivating primary source remains [Wood, arXiv:2501.16442v2,
section 4.1 and Appendix B](https://arxiv.org/pdf/2501.16442).
Its variable-interaction discussion does not supply the present analytic
prepared-sector theorem or a UV completion of this specified fixture.
The original fixture has actual `g`-coupled canonical clock matter and
free chi; no interacting-chi model is substituted for it.

The new proof is a direct coefficient-space argument from the literal
operator, not an application of an external effective-field-theory
theorem. Its elementary analytic ingredients agree with the
[NIST DLMF complex-variable calculus and power-series conventions,
section 1.9](https://dlmf.nist.gov/1.9). In particular the proof explicitly
derives the weighted-l1 inverse, derivative bounds and Taylor-tail norm;
it does not import an unquantified asymptotic expansion. The
Banach-valued Schwarz estimate is reduced to the scalar estimate using
bounded linear functionals and the dual norm. Neither DLMF nor the model
paper asserts this prepared-branch result.

The physical source normalization, conserved external TT stress and
fixed-slice endpoint maps are re-used only in their pinned S6.21 domains.
The improved transfer bound consequently retains `1<=K<=4`, whereas the
new analytic and center-coefficient arguments separately prove their
`0<=K<=4` scope, including the exact common mode at zero.

The action/sign dictionary is `+---`, `R_B=-6(DH+2H^2)` and
Einstein action `-M^2 R_B/2`; no photon-track opposite-curvature convention
is imported. The positive prepared kinetic factor belongs to a selected
source-free two-dimensional physical solution space at fixed K. It is not
a generally sourced spatially local EFT, a cone calculation, an
arbitrary-state stability theorem or a UV verdict. The source-dependent
and complex-delta proof assumptions are written explicitly in the
formulation and proof, not supplied silently by any citation.
