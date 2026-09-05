# Constant-Weyl candidate: finite-window tensor response

This new candidate adds cC integral C² in the unchanged physical matter
frame. Its old FLRW background is exact. The certificate controls the
chosen first-order, two-data tensor representative on a specified finite
hard-momentum window; it is not full candidate or quantum-causality control.

See the [contract](FORMULATION.md), [proof](notes/tensor-response.md) and
[certificate](certificates/cd-weyl-tensor.json).

## Replay

```bash
export PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src
.venv/bin/python -m p8_m1_weyl.verify --check
.venv/bin/python -m pytest -q problems/P8/s5/matter/physical/control/tree/loops/tensor/tests
```

The verifier is read-only and replays the pinned S5.9 lineage. Omit `--check`
to print the report. No finite matching coefficient, quantum vacuum,
fourth-order branch, scalar stability or UV interpretation is supplied by
this tensor comparison.
