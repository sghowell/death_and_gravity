# CD/M1 hard-channel tree control

This extension bounds the fixed matter-backed bounce's **cubic and quartic
finite-time tree transitions**, including matter mixing, on an explicit
high-frequency nonexceptional momentum domain. See the
[operational contract](FORMULATION.md) and [proof](notes/tree-bound.md).

The estimate uses exact coupled free modes and the full phase-space
canonical inverse. A positive coefficient recurrence bounds all retained
vertices and all four scalar/tensor internal channels without expanding
every quartic momentum kernel. The deterministic certificate computes a
single sufficient M*tau=10^324 for block norm at most 1/1000 at every center
time. That very loose scale is **not an optimized physical cutoff** or a
realistic duration proposal.

## Replay

From the repository root:

```bash
export PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src
.venv/bin/python -m p8_m1_tree.verify --check
.venv/bin/python -m pytest -q problems/P8/s5/matter/physical/control/tree/tests
```

The verifier is read-only; omit `--check` to print the reproducible report.
It replays and pins S5.7.CD, which checks the prior physical-reduction chain.

The scoped D-only and M1 tree-control gates now have parallel constructions.
Neither establishes higher-order/loop control, technical naturalness,
inclusive or forward scattering, a global vacuum, backreaction stability,
or UV admissibility. Those are the next distinct research obligations,
not consequences of small finite-order transition blocks.
