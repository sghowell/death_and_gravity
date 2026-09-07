# Convergent prepared tensor sector

`P8-S6.23.PREPARED` constructs a convergent analytic two-mode sector of the
unchanged full VARIABLE tensor system. It proves a positive **restricted**
physical kinetic normalization and a controlled source-free center
coefficient `1+(14/33)delta`, with error `<=48001delta^2`.

Fixed-slice regular-light error improves to an explicit `8600delta` upper
bound on the inherited spatial band. A genuinely conserved g-only source
with fixed temporal support can prepare the analytic sector; its required
delta-dependent retuning and finite cutoff norm are explicit.
This is not a general sourced EFT, a cone/causality calculation, or an
original-P8 matching/UV verdict.

See [the exact scope](FORMULATION.md), [the complete coefficient-space and
physical proof](notes/proof.md), and [source/sign boundaries](notes/sources.md).
The verifier pins and replays S6.21 and its full lineage without changing
any ancestor bytes. Independent Fraction arithmetic and a separately
authored full-action audit check the primary derivations.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/variable/response/prepared/tests
.venv/bin/ruff check problems/P8/s6/matching/variable/response/prepared
.venv/bin/python - --check <<'PY'
from pathlib import Path
import runpy
import sys
for source in Path('problems/P8').rglob('src'):
    sys.path.insert(0, str(source))
runpy.run_module('p8_variable_prepared.verify', run_name='__main__')
PY
```

Without `--check`, the verifier emits a candidate report to stdout. It
never writes or regenerates a certificate. Promotion requires independent
review; subsequent changes belong in a separately reviewed descendant.
