# Full spatial/projective BRST and fixed-source Ward bridge

This checkpoint continues the unchanged S259 parent. It supplies the
complete local classical spatial/projective ghost algebra, compensating
trace-gauge transformation, finite seven-ghost matrix and conditional
fixed-source mean/Hessian transport. It does not construct the full
interacting quantum state, regulator or mean.

Read [the formulation](FORMULATION.md), then
[algebra and canonical action](notes/algebra.md),
[projective compensation](notes/projective.md),
[translation and projection boundaries](notes/zero-modes.md),
[conditional Ward identities](notes/ward.md),
[fixed-state and source boundaries](notes/state.md),
[primary-source and original-problem scope](notes/scope.md), and
[immutable validation](notes/validation.md).

With all project source paths available, replay read-only using:

    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_spatial_brst.verify
    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_260

Only full regression uses the audited exact-GCD adapter:

    PYTHONPATH=problems/P8/s6/continuation/s6_219/tests \
      .venv/bin/python scripts/p8_replay.py full

No earlier scientific file or report is modified.
