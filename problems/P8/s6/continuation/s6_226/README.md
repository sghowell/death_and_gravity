# S6.226 — weighted original-channel inverse with bounded corrections

Original P8 remains OPEN. The original trace and pole-plus-cut shear factors now have a common lower bound on the full complex exterior. It gives a uniform exponentially weighted L2 inverse for a specified class of bounded time-dependent channel corrections, including their ordered composition with the actual S225 curvature coordinates.

This is not yet a bound on the actual unmatched curved correction. Removing the time weight can be extremely costly; a separate bounded correction with a growing reference pole demonstrates why weighted existence is not stability.

See [FORMULATION.md](FORMULATION.md) and the seven proof notes. Replay:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_226
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_weighted_channel_resolvent.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Native/direct/ordinary/CLI use original SymPy. Only full regression uses the audited exact-GCD adapter. All scientific, proof, test and native report bytes are frozen; post-freeze outcomes are recorded in the external audit.
