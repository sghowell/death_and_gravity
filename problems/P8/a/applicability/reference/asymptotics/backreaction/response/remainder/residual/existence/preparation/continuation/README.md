# A.13: an exact resummed causal block, not a long-time SEE solution

The full A.11 equation can be rearranged to treat its stiff Einstein term
and causal logarithm together. The resulting **frozen linear inverse**
has one growing pole, two damped poles and a cut. Its norm is not uniformly
small on longer intervals: for the named calibration it exceeds 1000 at
dimensionless duration 10^-5 and 10^28 at 10^-4. These intervals belong to
the comparator, not newly established durations of the actual solution.

The analysis retains every rolling/state term explicitly and proves the
exact growing-response cancellation condition for prescribed forcing.
It does not show that the actual nonlinear solution excites that pole,
remove it by a changed contour, or alter the physical prescription.

Read [FORMULATION.md](FORMULATION.md), the
[full-map dictionary](notes/decomposition.md),
[pole/inverse proof](notes/resolvent.md), and
[falsifiable next estimates](notes/next-step.md).

## Replay from the repository root

```sh
.venv/bin/pytest -q problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/tests
```

```sh
.venv/bin/python - <<'PY'
import runpy, sys
from pathlib import Path
root=Path("problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation")
sys.path[:0]=[str(root.parents[2]/"qsei"/"src"),*[str(p/"src") for p in (root,*root.parents)]]
sys.argv=["p8a_continuation.verify","--check"]
runpy.run_module("p8a_continuation.verify",run_name="__main__")
PY
```

The verifier checks the A.11 report and replays its complete lineage,
exact full-map/sign identities, a separate Fraction arithmetic engine,
source hashes and scope controls. It never writes its stored report.
The contour, pole-count, positivity and functional-analysis arguments are
written proofs, not a formalized global PDE theorem. All older sources
remain unchanged. P8(a) and P8 remain open.
