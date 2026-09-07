# Constant star-HR interactions: a regular flat no-bounce theorem

S6.16.STAR excludes physical contraction-to-expansion for the stated
constant pairwise star action, positive leaf Einstein coefficients, an
optional nonnegative central Einstein term, and a sole central NEC source.
It covers arbitrary interaction signs and Bianchi zero sets. It does not
assume subluminal cones or a nonsingular auxiliary TT Hessian.

The necessary exception is a disconnected nondynamical center: `G_u=0`
and all links endpoint-only. An actual algebraic solution also demonstrates
why algebraic leaves must not be inserted into an all-leg dynamic K.

Read [FORMULATION.md](FORMULATION.md), the
[proof](notes/proof.md), and the [primary-source audit](notes/sources.md).
The [certificate](certificates/star-no-bounce.json) is read-only; it pins
S6.14/S6.13 and the adopted S6 contract and hashes this subtree's direct
sources, tests and notes. It does not depend on the concurrent S6.15 child.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/star/tests
.venv/bin/python -c 'import pathlib,runpy,sys; sys.path[:0]=[str(p) for p in pathlib.Path("problems/P8").rglob("src")]; runpy.run_module("p8_star.verify",run_name="__main__")' --check
.venv/bin/ruff check problems/P8/s6/matching/star
```

Omitting `--check` emits a newly calculated candidate to stdout and writes
nothing. The verifier compares the complete report, including source hashes
and claim boundaries, not just selected output coefficients.

Core APIs: `potential.derive/evaluate`, `background.derive_two_leaf`,
`background.dynamic_reconstruction`, `branches.classify`, and the actual
`branches.algebraic_control/endpoint_control`. Numeric entry points accept
finite exact rational inputs (including rational strings), reject binary
floats and nonfinite values, and check their named positive domains.
`dynamic_reconstruction` is conditional on the supplied legs actually being
dynamic; it neither solves the complete equations nor chooses branches.
`compact_comparison_bound` accepts justified logarithmic-rate suprema,
not an unverified solution. The global no-crossing inference is a written
compact positive-part proof, not a numerical test of arbitrary functions.

This excludes a named family of candidate common parents. It supplies no
perturbative-health, finite-error matching or UV-completion verdict for the
original C/D rows, and does not close S6 or P8(b).
