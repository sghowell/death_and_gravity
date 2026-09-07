# Regular determinant-sum interactions: an actual-metric no-bounce theorem

S6.18.DETERMINANT treats the rank-one multi-vielbein determinant interaction,
with matter on actual Einstein metrics and separate physical NEC sources.
On the aligned positive flat-FLRW branch with both summed scale and lapse
nonzero, a positive weighted Einstein identity excludes contraction-to-
expansion. No high-frequency cone or scattering assumption is used.

Read the [formulation](FORMULATION.md), [full proof](notes/proof.md), and
[source audit](notes/sources.md). The
[certificate](certificates/determinant-no-bounce.json) pins S6.16, S6.13 and
the adopted S6 contract, not a mutable TREE checkpoint. It hashes the direct
new sources, tests and notes only.

The actual mixed-sign rolling/de Sitter fixtures verify background
equations, not perturbative health. A singular time-sum Bianchi control is
explicitly not a solution and prevents an invalid extension of the proof.
Moving matter to the auxiliary/summed metric would change the model.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/determinant/tests
.venv/bin/python -c 'import pathlib,runpy,sys; sys.path[:0]=[str(p) for p in pathlib.Path("problems/P8").rglob("src")]; runpy.run_module("p8_determinant.verify",run_name="__main__")' --check
.venv/bin/ruff check problems/P8/s6/matching/determinant
```

Without `--check`, the CLI prints a fresh candidate and writes nothing.
Whole-report equality checks all hashes and scope boundaries. Numeric
interfaces accept finite exact rationals and reject floats/nonfinite or
wrong-domain values. `potential.evaluate` reports singular sums without a
health verdict. `background.reconstruct` is explicitly conditional on an
already justified regular active component and physical-clock normalization;
it does not itself certify a solution or a sum's invertibility.

This excludes one specified parent family on a specified regular branch.
It does not close general matching, the original C/D rows, or P8(b).
