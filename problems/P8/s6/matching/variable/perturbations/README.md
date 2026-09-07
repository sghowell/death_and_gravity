# Regular scalar constraints and vector dynamics

New bounded claim `P8-S6.22.CONSTRAINTS`, using only the frozen S6.20
variable-beta background and its existing matching contract.

The result is a full quadratic scalar constraint reduction through the
bounce, an exact physical-observable dictionary, and a complete vector
action. It deliberately does **not** certify scalar stability: explicit
time/momentum-map controls expose why the apparent frozen scalar symbols
are insufficient.

- [Specified domain and exclusions](FORMULATION.md)
- [Written proof and normalization dictionary](notes/proof.md)
- [Primary-source applicability audit](notes/sources.md)
- [Read-only verifier](src/p8_variable_constraints/verify.py)

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/variable/perturbations/tests
```

The standalone CLI is `python -m p8_variable_constraints.verify --check`
with this and the pinned ancestor `src` directories on `PYTHONPATH`.
It does not write or repair the report. The report is
`certificates/variable-constraints.json`; it must be generated only after
independent review. The child does not edit or certify the separate
tensor-response or global variable-lapse work. Original P8 remains open.
