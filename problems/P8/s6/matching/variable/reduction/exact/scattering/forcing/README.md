# S6.36: one fixed pulse leaves a nonconvergent own-f response

Starting with zero hidden tensor data, one fixed smooth prescribed metric
pulse loads a provably nonzero state. The exact own-f equation carries it
through the central region with a phase-dependent original-Q readout:
its limsup/liminf gap is greater than eta/200 as delta tends to zero.
The physical duration and source shape are fixed.

This is a linear causal response on the unchanged off-shell physical
probe, not a conserved matter-source, low-frequency EFT or UV theorem.
It keeps both finite comparison errors; no convergent central remainder
or actual reference subsequential limit is assumed.

See [the statement](FORMULATION.md), [the full proof](notes/proof.md),
[the independent audit](notes/audit.md), and
[the source/domain interfaces](notes/interfaces.md). The source-pinned
report recursively rebuilds S6.35 and its unchanged ancestry.

## Ordinary replay

```sh
PYTHONHASHSEED=0 .venv/bin/python -c 'from sympy.core.random import seed; seed(0); import pytest; raise SystemExit(pytest.main(["problems/P8/s6/matching/variable/reduction/exact/scattering/forcing/tests", "-q", "-p", "no:faulthandler"]))'
```

With local P8 src directories on the Python path, the separate CLI is
`python -m p8_own_forced.verify --check`. Without --check it prints the
report and writes no files. Ordinary certificate tests and CLI use no
broad-regression exact-GCD adapter. Both seeds are zero and the diagnosed
host timer plugin is disabled.

Original P8 remains open. The scoped P8(a) photon objective is unchanged.
