# Full physical tensor response on fixed slices

`P8-S6.21.RESPONSE` extends the frozen VARIABLE tensor calculation from a
shrinking inner limit to the full four-component propagator between the
fixed physical slices `T=+-tau/100`, uniformly for `1<=tau*kcom<=2`.
All physical normalization, moving-weight, source and phase terms remain.

The adapted physical Cauchy error is at most `40000 delta^(1/3)` for
`0<delta<=10^-9`; a nontrivial example is `<=1/250` for `delta<=10^-21`.
An explicit conserved g-only source construction realizes either incoming
relative-mode data or prepared regular-light data, with its exact source
norm retained in the finite-parameter bound. Neither fixed duration nor a
fixed spatial momentum band is a certified low temporal-frequency band.
This is not a light-only EFT or original-P8 UV verdict.

See [FORMULATION.md](FORMULATION.md), [the full norm proof](notes/proof.md),
and [source/sign boundaries](notes/sources.md). The certificate pins and
replays VARIABLE and its adopted-contract lineage without modifying them.
The continuous coefficient enclosures use independent Fraction intervals
and bivariate Taylor jets; an independently authored physical-action audit
is included. Exact arithmetic and the written ODE argument are the proof
boundary, not a proof-assistant formalization.

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/variable/response/tests
.venv/bin/ruff check problems/P8/s6/matching/variable/response
.venv/bin/python - --check <<'PY'
from pathlib import Path
import runpy
import sys
for source in Path('problems/P8').rglob('src'):
    sys.path.insert(0, str(source))
runpy.run_module('p8_variable_response.verify', run_name='__main__')
PY
```

The verifier is read-only. Without `--check` it emits a candidate report to
stdout; it never creates or overwrites the pinned certificate. New work
belongs in a separately reviewed descendant, not in frozen source files.
