# S6.243 — complete heavy scalar spatial UV input and finite difference

The unchanged QG2 minimally coupled heavy scalar now has its COMPLETE local spatial UV input. Every ordered trace/gradient product and source endpoint is retained. The full six-invariant logarithmic symbol matches the scalar covariant pole, and the entire dimension-dependent measure and counteraction give the actual same-scheme finite UV difference.

All eighteen finite coefficients obey a normalized bound below10^-600 on detector L2 times source Z24, with two time and four spatial derivatives, on |t|<=1/2 and all spatial transfers. This is not the full UV-subtracted comparison response, homogeneous anchor, heavy ADM/clock response or a quantum inverse. Original P8 remains OPEN.

See [the formulation](FORMULATION.md), [full vertex](notes/vertex.md), [ordered endpoints](notes/ordering.md), [dimensional geometry](notes/dimension.md), [full finite prescription](notes/finite.md), [complete local bound](notes/bounds.md), [scope](notes/scope.md) and [validation](notes/validation.md).

Read-only checks:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_243
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_heavy_spatial_uv.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Native/direct/ordinary/CLI retain original SymPy. Only full regression uses the audited exact-GCD adapter. Frozen scientific inputs are unchanged. The written proofs are not FORMALIZED.
