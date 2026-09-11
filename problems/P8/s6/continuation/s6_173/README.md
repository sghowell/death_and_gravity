# P8 S6.173: leading actual CD clock sources

The full analytic target's first clock jets give an exact
all-time source decomposition, with original M1 matter retained.
The leading curvature, A3 and scalar contributions are nonzero;
A4,A5 remain necessary to the constraint block despite their zero
first background source. The old vacuum-unit action norm remains
outside the large-gradient clock domain.

See FORMULATION.md and notes/ for the exact scope. This is not a
propagating parent or quantum bounce; original V/G/B and P8 are OPEN.

After publication, from the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_173
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_target_clock_null.verify

Native/direct/ordinary/CLI keep original SymPy. Only the separately
audited full regression uses its exact GCD adapter.
