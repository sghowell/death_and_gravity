# S6.196: separated-support spatial continuum response

The actual Gaussian Proca current now has a full-continuum weak
first-order response on strictly time-ordered smooth compact
spatial shear supports, with a quantitative regulator-removal error.
This does not extend the retarded kernel to overlapping times or
the diagonal. Original P8 remains OPEN.

Read FORMULATION.md and the seven notes.

## Read-only replay

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_196
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_separated_response.verify
    .venv/bin/python -u scripts/p8_replay.py full

All eighteen scientific sources and the native report are hash-pinned.
Native, direct science, ordinary and CLI use original SymPy.
Only full regression uses the audited exact GCD adapter.
Continuous Hilbert-space and support arguments are not FORMALIZED.
