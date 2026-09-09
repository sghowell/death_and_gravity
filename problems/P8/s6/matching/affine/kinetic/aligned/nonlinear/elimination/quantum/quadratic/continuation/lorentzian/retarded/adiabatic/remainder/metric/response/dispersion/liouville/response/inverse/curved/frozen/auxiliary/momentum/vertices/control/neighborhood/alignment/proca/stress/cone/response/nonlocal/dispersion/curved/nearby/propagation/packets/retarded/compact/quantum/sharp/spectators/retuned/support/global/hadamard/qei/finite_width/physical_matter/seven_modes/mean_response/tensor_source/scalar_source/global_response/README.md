# Global coupled scalar and joint seven-mode relative response

The scalar infinite-tail estimate is obtained from the actual
S6.106 kernels. No minimally coupled scalar source is substituted.
The damping-weighted phase and compensated trace expose
integrable actual source kernels without altering physical means.

Read the [formulation](FORMULATION.md), [phase](notes/phase.md),
[compensation](notes/compensation.md), [global bounds](notes/bounds.md),
[physical observable](notes/observable.md), [joint state and response](notes/joint.md),
[geometry](notes/geometry.md) and [scope](notes/scope.md).

With the repository environment and all P8 src directories
on sys.path, native read-only replay is

    python -m p8_scalar_global_response.verify --check

The separately audited complete regression is

    .venv/bin/python -u scripts/p8_snapshot_regression.py

Only the latter full-regression runner uses its exact GCD
adapter. Native science, report generation, ordinary replay
and the CLI do not monkeypatch scientific SymPy.
Every own source/proof/test byte and every report field is
pinned and independently checked.

The [report](certificates/global-seven-mode-relative-mean-control.json)
reconstructs complete rational kernels before computing exact
monomial fingerprints. It contains explicit coefficientwise
envelopes and exact rational bounds, not momentum/time samples.
Written analytic proofs are not proof-assistant formalized.
Original P8 remains open.
