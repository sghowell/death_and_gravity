# P8 S6.222: quantum-forced scalar interface

The complete forced constraints contain quantum-force corrections in both lapse and shift. Their exact phase/metric maps give an ordered Schur equation on a stated graph domain. It is not a bounded quantum inverse or P8 closure.

Read [FORMULATION](FORMULATION.md), [forces](notes/forces.md), [feedback](notes/feedback.md) and [domains](notes/domains.md). The finite-dimensional controls establish algebra/order checks, not the continuum inverse.

Ordinary: `.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_222`.
CLI: `.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_quantum_forced_constraints.verify`.
Full: `PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full`.
Only full regression uses the audited exact-GCD adapter. The helper path leaves the frozen S219 test unchanged.
