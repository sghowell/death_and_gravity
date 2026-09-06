# P8(a) A.12: a QSEI on the actual local semiclassical solution

The new exact A.11 solution satisfies an absolute all-Hadamard,
all-H2-sampler QSEI with sufficient coefficient `2 hbar/(16pi^2)`
on its source-free half-slab. A new two-frequency argument needs only
the known potential and first-derivative bounds. Positivity of the
reference EED is derived from the actual SEE.

The available interval still cannot meet the cosmological comoving
focusing test. This is not a singularity theorem or completion of P8.
Read the [formulation](FORMULATION.md), [proof](notes/proof.md) and
[quantitative focusing boundary](notes/focusing.md).

From the repository root:

```sh
.venv/bin/python -m pytest problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei/tests -q
.venv/bin/python - <<'PY'
import runpy
import sys
from pathlib import Path

root = Path("problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei")
sys.path[:0] = [str(root.parents[2]/"qsei"/"src"),
               *[str(path/"src") for path in (root, *root.parents)]]
sys.argv = ["p8a_see_qsei.verify", "--check"]
runpy.run_module("p8a_see_qsei.verify", run_name="__main__")
PY
```

The CLI is read-only. Without `--check` it prints a candidate report,
never overwriting stored evidence. The pinned certificate is
`certificates/see-qsei.json`. No download or external service is needed
for replay.
