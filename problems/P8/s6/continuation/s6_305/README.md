# P8 S6.305: fixed spectator metric insertions and known-reference rate

This immutable successor assembles the entire source-fixed H/Proca Gaussian
metric kernels in all crossed physical four-Phi channels and adds the
original M1 logarithmic part once. It retains the light-Phi loop already in
S302 without double counting. It derives both optical normalizations,
the same-axis finite matching map, a uniform signed known-reference rate
and its forward Newton limit. Original V/G/B/P8 remain OPEN.

Read [the contract](FORMULATION.md), [source ownership](notes/source.md),
[the insertion](notes/insertion.md), [cuts and sheets](notes/cuts.md),
[uniform bounds](notes/bounds.md), [matching](notes/matching.md),
[scope](notes/scope.md) and [validation](notes/validation.md).

The native20-field report is
`certificates/polynomial-vacuum-affine-spectator-gravity-insertion.json`.
Own replay: `.venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_305`.
CLI replay: `.venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_spectator_gravity_insertion.verify`.
Full acceptance also requires the captured complete P8 regression.
