# Constant-Weyl CD/M1 coupled scalar response

Read [FORMULATION.md](FORMULATION.md) and the
[written proof](notes/scalar-response.md) before interpreting the certificate.
This new subtree continues the same constant-Weyl candidate after S5.10's
tensor result; it does not modify the frozen CD/M1 action.

From the repository root:

```sh
.venv/bin/python -m pytest -q problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/tests
```

For replay, put this subtree and its ancestor packages on `PYTHONPATH`:

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/matter/src:problems/P8/s5/matter/physical/src:problems/P8/s5/matter/physical/control/src:problems/P8/s5/matter/physical/control/tree/src:problems/P8/s5/matter/physical/control/tree/loops/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/src:problems/P8/s5/matter/physical/control/tree/loops/tensor/scalars/src .venv/bin/python -m p8_m1_weyl_scalar.verify --check
```

The read-only verifier replays the pinned lineage, identities, negative
controls and exact finite-band margins. It hashes this subtree's Python
sources, tests and written proof. Certificate generation prints JSON; `--check`
does not write any file. No higher-derivative branch or quantum causality
conclusion follows from a successful replay.
