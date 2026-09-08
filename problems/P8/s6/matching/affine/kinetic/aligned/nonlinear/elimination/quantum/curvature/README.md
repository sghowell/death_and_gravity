# Local Proca curvature checkpoint

[Scope](FORMULATION.md), [local/physical-geometry proof](notes/local.md),
[independent spectral control](notes/spectral.md).

S6.49 computes the retained vector's on-clock local curvature pole
and global physical coefficient bounds. A separate bounded S4
spectral calculation checks its normalization and excluded-mode
bookkeeping. Original P8 remains open: this is not a finite curved
or in-in quantum error estimate, nor a quantum-corrected bounce.

The ordinary tests are in `tests/`; their conftest discovers all
P8 source roots without changing ancestor files. Use the seeded
Python/pytest recipe in the current audit, with `-p no:faulthandler`.
After exposing those source roots, `python -m p8_vector_curvature.verify
--check` performs read-only replay. The CLI never writes a report.
