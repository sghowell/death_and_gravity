# P8(a) A.10: actual response continuity and conditional inverse tools

The actual nonlinear massless mode functional is now quantitatively
Lipschitz in the uniform norm of the potential derivative. The proof
needs no frequency gap or higher potential jets and preserves a shared
nonzero preparation history. A separate causal logarithmic inverse has
an explicit positive pole/cut formula and rational norm bounds.

These are ingredients, **not an actual semiclassical existence theorem**.
Read the [formulation](FORMULATION.md), [mode proof](notes/mode-response.md),
[inverse proof](notes/inverse.md) and [remaining obligations](notes/see-gap.md).

From the repository root:

```sh
.venv/bin/python -m pytest problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/tests -q
.venv/bin/python - <<'PY'
import runpy
import sys
from pathlib import Path

root = Path("problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence")
sys.path[:0] = [str(root.parent/"qsei"/"src"),
               *[str(path/"src") for path in (root, *root.parents)]]
sys.argv = ["p8a_existence.verify", "--check"]
runpy.run_module("p8a_existence.verify", run_name="__main__")
PY
```

The second command adds this subtree, its QSEI sibling and ancestor
packages to the import path, as the test `conftest.py` does.
The certificate CLI is read-only; omitting `--check`
prints a fresh candidate report and never overwrites stored evidence.
No downloads, paid compute or external services are needed for replay.
