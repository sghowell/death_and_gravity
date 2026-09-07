# Source-pinned local vacuum prerequisite

This child tests the actual constant-clock vacuum prerequisite of the same
variable-beta action used for S6.20's local CD-shaped background. It does
not replace that action by its value at a moving background metric.

The independently reviewed continuous proof finds no positive proportional Minkowski
vacuum anywhere in the specified local clock interval, for either admitted
action-parameter branch. The independent vacuum ratio is unconstrained a
priori. Its f-flat root leaves a strictly negative g-flat residual; the
literal clock equation is also derived and retained.

See [FORMULATION.md](FORMULATION.md), the [proof](notes/proof.md), and
[source and scope audit](notes/sources.md). The
[certificate](certificates/local-vacuum-obstruction.json) pins all local
sources, including the separately authored audit, and recursively replays
the actual S6.20 action. From the repository root, a fresh ordinary exact
replay is:

```sh
PYTHONHASHSEED=0 .venv/bin/python - <<'PY'
import runpy
import sys
from pathlib import Path
from sympy.core.random import seed
seed(0)
for source in Path("problems/P8").rglob("src"):
    sys.path.insert(0, str(source.resolve()))
sys.argv = ["p8_variable_vacuum.verify", "--check"]
runpy.run_module("p8_variable_vacuum.verify", run_name="__main__")
PY
```

The Python hash seed does not seed SymPy's separate factorization RNG,
so both are set explicitly. No arithmetic adapter or faulthandler timer
is used. For pytest, also disable its faulthandler plugin with
`-p no:faulthandler`; the host's diagnostic watchdog can otherwise hang
independently of the proof computation.

Smooth off-interval extensions and controlled matching remain separate
research tasks. No original-P8-wide, all-parent or UV exclusion is asserted.
