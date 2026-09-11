# P8 S6.174: regular classical affine/vector parent

CD-REG-AFFINE-ISO gives a separately specified classical common-domain
construction on the full analytic CD target. Exact auxiliary reduction
retains that target and a sourced Proca field. The source vanishes
through first order on the entire clock. Both the vacuum and clock
have the stated healthy quadratic blocks; local nonlinear seven-mode
constraint neighborhoods exist at each finite clock time.

Quantum matching, cutoff, controlled background/response and original
V/G/B and P8 remain OPEN. See FORMULATION.md and notes/ for the
precise action, boundaries and scope.

Read-only replay from the repository root after publication:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_174
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_analytic_affine_parent.verify

Native/direct/ordinary/CLI use original SymPy. Only the separately
audited full regression runner uses the exact GCD adapter.
