# S6.26.FORCED: causal actual-source tensor response

Reviewed mathematical gate with a source-pinned, read-only certificate.

For the same VARIABLE parent, zero physical initial data and arbitrary
bounded conserved external `g`-TT forcing, this gate gives an explicit
retarded-memory approximation to the **physical g response**:

`||g-gapp||infinity <=1000*r^4*||sigma||infinity`.

The domain is `0<r<=1/100`, `0<delta<=min(10^-9,r^2/100)`, `0<=K<=4` on
`u in [-r,r]`. The actual source and delta-dependent physical map are retained.
At `r=1/100` the error is at most `10^-5` times the source supremum, uniformly
in delta. It is not a relative error or delta-convergence statement.

The approximation retains the exact universal heavy retarded kernel; it is
not a local EFT. The source-free analytic sector has a different projected
source factor. A separate Poincare diagnostic limits only a specifically
defined RMS-time-scale / instantaneous-stiffness-proxy hierarchy for compact
pulses. Neither diagnostic supplies a rolling gap, cutoff or UV verdict.

Read [FORMULATION.md](FORMULATION.md), the [full proof](notes/proof.md), and
the [primary source/scope audit](notes/sources.md). Public modules expose the
literal operator and source, exact kernel, rational continuous envelopes,
finite error API, and independent Fraction reconstruction. No frozen
ancestor has been edited.

From the repository root:

```sh
.venv/bin/pytest -p no:faulthandler -q problems/P8/s6/matching/variable/response/forced/tests
.venv/bin/ruff check problems/P8/s6/matching/variable/response/forced
PYTHONHASHSEED=0 .venv/bin/python - <<'PY'
from pathlib import Path
import runpy
import sys
from sympy.core.random import seed
seed(0)
for source in Path('problems/P8').rglob('src'):
    sys.path.insert(0, str(source))
sys.argv = ['p8_variable_forced.verify', '--check']
runpy.run_module('p8_variable_forced.verify', run_name='__main__')
PY
```

The final command is read-only: it compares the frozen certificate with a
fresh exact reconstruction and never regenerates or overwrites the report.
Both Python's hash seed and SymPy's polynomial-algorithm RNG are fixed for
reproducible recursive replay. The pytest faulthandler plugin is disabled
because this host's watchdog cancellation can spin; no faulthandler timer is
started. No cached prior result or arithmetic adapter replaces replay.
The verifier pins S6.23 and recursively rebuilds its frozen S6.21/S6.20
lineage, checks the original S6 contract, and compares the report exactly.

This gate does not close original S6 or P8.
