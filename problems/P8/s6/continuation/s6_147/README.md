# Finite MS quadratic primitive mass references

The massless-exchange anchors are computed from full dimensional
mass tensors and proper counterterms. A fixed-scale vacuum
derivative provides an independent normalization. The actual scalar
mass is restored by a regulated integral over its mass squared;
the proper quartic counterterm and all finite pole products remain.

Run the ordinary and command-line checks from the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_147
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_fermion_ms_mass.verify

Only the full P8 regression uses its separately audited GCD adapter.
Native, ordinary, CLI and direct science use unmodified SymPy.
See FORMULATION.md and the six written proof notes. Exact tests are
not formalization or independent peer review. Original P8 remains OPEN.
