# S6.30: the own-f action does not yet match original CD

The unchanged variable-coupling two-metric parent has a well-defined
**formal** stationary-f derivative expansion. Through four total
derivatives it gives a quartic-Horndeski scalar-tensor part plus a
nonzero curvature-square operator, with the physical g and free chi
left unchanged. It is not the original CD/M1 covariant action.

Three independent tests distinguish this result from a background match:

- The literal physical tensor equation of the retained action has a
  fourth-order term. Every original quadratic-DHOST action is second
  order on the same off-shell tensor test.
- In the specified scalar-tensor normal form, the clock-invariant
  quantity X(A1-2F2_X)/G_T is0 at the bounce; CD requires1 throughout
  the open clock tube. Exact matching needs a normalized F2_X remainder
  of magnitude1/2 even after a favorable leading Planck normalization.
- The actual background lapse defect tends to zero at the center, but
  its fourth physical-clock derivative tends to10752. The first metric
  correction has a nonvanishing second-derivative limit. Small values
  do not imply a uniformly small derivative expansion.

Removing the curvature-square term by a metric redefinition generates
new chi-curvature interactions and, after a further specified
leading-equation rewrite, chi self-interactions. Their sources and
contacts cannot be discarded while claiming the same free-M1 frame.

These are exact retained-action and necessary-remainder statements.
There is no proved convergence or small omitted remainder, rolling gap,
physical cutoff, parent ghost, UV verdict or whole-row exclusion.
Higher-order/resummed descriptions need their own controlled dictionary.
The word formal refers to derivative grading, not proof-assistant
formalization. Original P8 remains open.

See [the precise statement](FORMULATION.md),
[stationary action and background jets](notes/stationary.md),
[the clock/matter dictionary](notes/dictionary.md),
[the independent tensor test](notes/tensor-obstruction.md), and
[source and verification boundaries](notes/sources.md).

## Read-only reproduction

From the repository root, using the installed virtual environment:

```sh
PYTHONHASHSEED=0 .venv/bin/python -c '
from sympy.core.random import seed
seed(0)
import pytest
raise SystemExit(pytest.main([
    "problems/P8/s6/matching/variable/reduction/tests", "-q",
    "-p", "no:faulthandler"]))
'
```

For a separate fresh ordinary certificate replay:

```sh
PYTHONHASHSEED=0 .venv/bin/python -c '
import sys, runpy
from pathlib import Path
from sympy.core.random import seed
seed(0)
for path in Path("problems/P8").rglob("src"):
    sys.path.insert(0, str(path))
sys.argv=["p8_variable_reduction.verify", "--check"]
runpy.run_module("p8_variable_reduction.verify", run_name="__main__")
'
```

The replay rebuilds the frozen parent recursively and verifies the
original CD witness and its source hashes. It does not trust an old
passed flag. The separately seeded factor RNG is deliberate; the host's
problematic faulthandler timer is disabled without removing assertions.
The report is [stationary-cd-matching.json](certificates/stationary-cd-matching.json).
