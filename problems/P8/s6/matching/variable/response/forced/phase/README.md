# S6.27.PHASE: fixed-source physical phase separation

Reviewed gate with a [source-pinned read-only certificate](certificates/fixed-source-phase-separation.json).

For every fixed `0<r<=1/1000`, fixed `0<=K<=4`, and the stated class of
delta-independent near-rectangular pulses, including compact smooth
pulses, the full physical endpoint response is bounded but has no
delta->0 limit. Its limsup minus liminf is greater than `3*r^2/20`.

The proof keeps the actual source/map and bounds the full coupled
remainder. It does not infer the result from an impulse diagonal or
sampled oscillations. It excludes this fixed-source limiting response,
not all local EFTs or a delta-dependent memory description.

Read [FORMULATION.md](FORMULATION.md), the [complete proof](notes/proof.md)
and the [primary-source audit](notes/sources.md). `core.separation_bound`
only calibrates a conditional theorem; numeric parameters cannot certify
an arbitrary source profile's functional hypotheses.

From the repository root:

```sh
.venv/bin/pytest -q -p no:faulthandler problems/P8/s6/matching/variable/response/forced/phase/tests
.venv/bin/ruff check problems/P8/s6/matching/variable/response/forced/phase
PYTHONHASHSEED=0 .venv/bin/python - <<'PY'
from pathlib import Path
import runpy
import sys
from sympy.core.random import seed
seed(0)
for path in Path('problems/P8').rglob('src'):
    sys.path.insert(0, str(path))
sys.argv = ['p8_forced_phase.verify', '--check']
runpy.run_module('p8_forced_phase.verify', run_name='__main__')
PY
```

The final command checks the reviewed frozen report; it rebuilds and
compares without writing. Both RNGs are fixed and no faulthandler timers
are used. No frozen ancestor or original P8 verdict is changed.
