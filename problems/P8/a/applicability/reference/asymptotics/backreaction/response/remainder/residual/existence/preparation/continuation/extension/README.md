# A.14: the same actual solution on a 10,000-times longer interval

The original preparation and state now have a proved smooth continuation
to dimensionless future length 10^-6. The source still turns off at
5*10^-11; its support, cutoff, past, physical prescription and radiation
normalization are unchanged. Uniqueness proves agreement with A.11 on
the original interval.

The full weighted fixed-point map includes every Einstein, source,
anomaly, local and actual mode-response term. Its contraction is below
3*10^-6. Separate pointwise barriers survive the exp(2) conversion; this
is an actual full SEE result, not merely a weighted linear estimate.

See [FORMULATION.md](FORMULATION.md), the
[complete contraction proof](notes/contraction.md), and
[smoothness/state/constraint proof](notes/regularity.md).
The extended-domain QSEI and focusing check remain separate tasks; P8
is not closed.

## Read-only replay from the repository root

```sh
.venv/bin/pytest -q problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/extension/tests
```

```sh
.venv/bin/python - <<'PY'
import runpy, sys
from pathlib import Path
root=Path("problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/extension")
sys.path[:0]=[str(root.parents[3]/"qsei"/"src"),*[str(p/"src") for p in (root,*root.parents)]]
sys.argv=["p8a_extension.verify","--check"]
runpy.run_module("p8a_extension.verify",run_name="__main__")
PY
```

The stored report hashes the new source/proof/tests, pins A.13 and A.11,
and replays their complete prior lineage plus independent Fraction
arithmetic. The CLI never writes a report. Functional-analytic and
Hadamard implications remain written proofs, not finite mode samples
or a formalized global PDE theorem.
