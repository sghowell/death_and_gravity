# S6.272: finite-regulator volume turnaround

For the same full-source, 48-pair finite regulator, this packet extends the
evaluated real interval to |u| <= 10^-130 and proves a state-independent
volume turnaround for both calibrated-coherent and Weyl orderings and both
original cutoffs. Original P8 remains OPEN.

The actual Weyl normalized volume operator differs from I by less than
10^-264 and is greater than I/2. Multiplying by the unchanged physical
reference volume (1+u^2)^6, both endpoint means exceed the center mean by
more than 5 x 10^-260. A global minimum therefore occurs in the interior;
every minimizer satisfies |u| < 10^-132.

This is a minimum of an explicitly regulated readout with EXTERNAL
homogeneous background functions. It is not a unique or strict bounce,
self-consistent quantum background evolution, an unlocalized interacting
mean, a continuum limit, matching, omitted-loop control, or UV completion.
No small longer-time state displacement or phase-space leakage is proved.

Read FORMULATION.md and all six notes for the complete written proof.
The certificate checks exact source identities, rational bounds, source
hashes and ancestry; it does not formally verify the analytic proof.
Independent full-metric and quantum fixtures are diagnostics, not a
simulation defining the full P8 dynamics.

Read-only replays:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_272
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_finite_volume_turnaround.verify
