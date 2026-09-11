# S6.195: spatial two-momentum Proca current

The unchanged conditional vector sector has explicit full spatial
metric vertices and contacts, a two-momentum Gaussian tangent,
and quantitative finite-band bounds uniform in external momentum.
This is not a full renormalized spatial continuum response or
a physical cutoff. Original P8 remains OPEN.

Read FORMULATION.md and the seven notes.

## Read-only replay

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_195
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_spatial_current.verify
    .venv/bin/python -u scripts/p8_replay.py full

All eighteen scientific source files and the native report are
hash-pinned. Native, direct science, ordinary and CLI use original
SymPy. Only the complete regression uses the audited exact GCD
adapter. Continuous arguments are not FORMALIZED.
