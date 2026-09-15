# P8 S6.275: canonical-boundary correction and new finite hybrid

Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

This successor corrects the S269 raw/prepared scalar momentum dictionary:
the earlier S220 integration-by-parts boundary was missing. At the bounce,

raw r = prepared pv + 3 sigma/10,
raw m = prepared ps + 3 v/10.

It changes momentum correlations of the original prepared Gaussian and
is NOT a removable global scalar phase. The frozen S269-S274 numerical
models/reports remain archived, but their same-original-physical-state
identification is explicitly qualified. None is silently rewritten.

After retaining the exact boundary shear and its time generator, the
same original prepared seed/source and regulator prescriptions give
FOUR NEW corrected finite classical-homogeneous/quantum-mode solutions.
These are not retrospectively identified with the archived trajectories.

On |u|<=1e-180, the corrected coupled model has a unique solution in the
declared domain for each of the two cutoffs and complete orderings.
Homogeneous deviation is<1e-390, volume relative error<1e-380, endpoint
gap>5e-360 and all minima lie in |u|<1e-188. The positive coherent
outside-core probability is<1e-6. Coupled cutoff(Y,phase-factored state,
volume) comparisons are<1e-580,1e-12,1e-514; ordering comparisons
are<1e-630,1e-65,1e-566.

These finite results do not prove sharp support, a unique/strict
minimum, homogeneous quantization, regulator removal, physical UV
matching, omitted-loop/Regge control or nonlinear global completion.

Read FORMULATION.md, then notes/boundary.md for the physical correction.
The other notes state complete source, operator, dynamics and scope
arguments. The exact audit has111 named checks,1493 scalar entries,
94 proof gates and13 controls including623 rejected inputs.
All9 original primitive and130 historical matching rows are retained;
H8A439 records the correction and new bounded result, with explicit
qualification of H8A433-H8A438.

Read-only validation from the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_275
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_canonical_boundary_corrected_hybrid.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native/direct/ordinary/CLI require original SymPy. Only the audited FULL
runner may use its exact-GCD adapter. A replay verifies frozen contents;
it cannot erase a separately proved physical erratum.
