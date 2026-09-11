# S6.194: complete fixed-prescription homogeneous current

This checkpoint gives a quantitative bound for the actual
Gaussian Proca homogeneous current, its first two shear-amplitude
derivatives, its complete finite-amplitude Taylor remainder and
the correctly normalized linear response at the isotropic clock.
It retains the same state and finite counterterms.

Read FORMULATION.md and the seven written proof notes. This is
not a full spatial/mixed response inverse or interacting quantum
background; original P8 is OPEN.

## Read-only replay

From the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_194
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_covariant_current.verify
    .venv/bin/python -u scripts/p8_replay.py full

The native report and all eighteen scientific source files are
hash-pinned. The ordinary and CLI replays use original SymPy;
only the complete captured regression uses the audited exact
GCD adapter. Written continuous arguments are not FORMALIZED.
