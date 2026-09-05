# CD/M1 isolated-matter one-loop audit

This extension computes the mandatory curvature-squared local counterterm
of the free real M1 scalar and proves strict finite-basis nonclosure in the
frozen physical matter frame. It also bounds the specified local logarithmic
source on CD; it does **not** bound the complete quantum correction.

See the [contract](FORMULATION.md), [derivation and source conventions](notes/one-loop.md)
and [certificate](certificates/cd-matter-loop.json). S6 permits the higher
operators as new EFT matching candidates. No old certificate is altered.

## Replay

From the repository root:

```bash
export PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src
.venv/bin/python -m p8_m1_loops.verify --check
.venv/bin/python -m pytest -q problems/P8/s5/matter/physical/control/tree/loops/tests
```

The verifier is read-only; omit `--check` to print the report. It pins and
replays S5.8 and its complete prior chain. Independent tests retain the
lapse until after variation and directly contract the linearized Weyl
tensor, checking two misleading background-only shortcuts.

Finite matching coefficients, massless nonlocal and state-dependent effects,
other loops, corrected causality/backreaction, higher-order errors and UV
matching remain open. The enormous sufficient S5.8 scale is not an optimal
cutoff or a bound on those missing effects.
