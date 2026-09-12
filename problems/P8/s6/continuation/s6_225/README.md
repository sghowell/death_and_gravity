# S6.225 — actual curved scalar geometry and ordered reference inverse

Original P8 remains OPEN. This checkpoint derives the actual curvature-coordinate matrix on the unchanged background, retains the entire original finite local scalar Hessian, and constructs a uniformly bounded inverse for a precisely specified curvature-adapted reference.

The reference uses the S224 flat trace and shear scalar factors. It is not an identification with the complete curved quantum state response. The explicit finite local remainder, curved mass/state/contact/tree/matter terms and full S222 auxiliary/clock interface are not discarded or claimed small.

Read [FORMULATION.md](FORMULATION.md) and the seven written proof notes. Replay:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_225
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_curved_scalar_reference.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Native, direct, ordinary and CLI runs retain original SymPy. Only full regression uses the audited exact-GCD adapter. The report, every scientific source and every proof/test input are frozen and hash-checked. Fresh post-freeze outcomes belong in the external audit, not in rewritten frozen sources.
