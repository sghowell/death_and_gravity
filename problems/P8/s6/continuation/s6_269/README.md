# P8 S6.269: evaluated finite nonlinear bounce phase domain

This continuation evaluates an explicit finite member of S268's
local coherent quantum regulator. At u=0, take the L=1 torus,
wavevectors +-10^64 e_i for all three axes, all eight channels
and the SAME original pure finite nonzero-mode seed.

The entire96-dimensional whitened phase ball of radius2*10^20
reconstructs full nonlinear spatial fields and twelve actual
auxiliary invariants below10^-260. The complete source,
cotangent/ghost corrections, curvature, density factors and
generated Fourier harmonics remain. The coherent core of
radius10^20 has strictly positive leakage below exp(-10^39).

The positive finite coherent physical volume obeys
||Q_V(Fext)-I||<10^-255. For two admissible core-agreeing
volume cutoffs, the same-seed vector and mean differences
are bounded by2*10^-255 exp(-10^39/2) and
2*10^-255 exp(-10^39), respectively.

This is a BOUNCE-SLICE regulated-observable result. It does
not evaluate a positive time interval, interaction/ordering
error, calibrated volume, original interacting mean,
physical cutoff or continuum limit. Original V/G/B/P8
remain OPEN. Completed scoped P8(a) is unchanged.

Read FORMULATION.md and the seven written-proof notes.
The read-only verifier rebuilds every frozen ancestor;
the certificate contains exact rational identities and
bounds plus scoped written proof gates, not a Lean proof.

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_269
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_quantitative_phase_domain.verify

Native/direct/ordinary/CLI use original SymPy. The separate
full regression alone uses the audited exact-GCD adapter.
The external publication assessment records observed
exit statuses and raw source/report validation.
