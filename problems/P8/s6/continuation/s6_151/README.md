# Wineglass MS interaction forest

Read [FORMULATION.md](FORMULATION.md) and the six proof notes.
The finite equal-mass constant is represented by a convergent
integral with an analytic enclosure, not by its diagnostic
decimal. Both the local reference and full nonlocal scale term
are included in the family bound.

The read-only report pins eighteen source/proof/test files and
the frozen S6.150 parent. From the repository:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_151
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_scalar_zero_reference.verify

Only full regression uses the audited exact GCD adapter.
Native, ordinary, CLI and direct science use unmodified SymPy.
This is one completed interaction-family conversion, not the
complete scalar/GY14 MS matching or original P8 closure.
