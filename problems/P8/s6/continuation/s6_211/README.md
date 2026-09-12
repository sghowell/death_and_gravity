# S6.211: actual spatial UV difference and fixed pole identity

The complete actual unit-W8 endpoint symbol is compared at transfer P and P0, retaining the fixed homogeneous anchor. The full one-leg contact cancels only in that spatial difference. Exact curved logarithmic coefficients agree with the original pole Hessian.

This does not complete finite dimensional matching, the full quantum response/inverse, background, remaining V/G/B gates or original P8.

See [formulation](FORMULATION.md), [jets](notes/jets.md), [spatial coefficients](notes/density.md), [contact](notes/contact.md), [pole boundary](notes/matching.md), [norms](notes/norms.md), [scope](notes/scope.md) and [verification](notes/validation.md).

Read-only frozen replay commands from the repository root:

    .venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_211
    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_spatial_matching_difference.verify
    .venv/bin/python -u scripts/p8_replay.py full

Native, direct, ordinary and CLI calculations retain original SymPy. The audited exact GCD adapter is used only for full regression.
