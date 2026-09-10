# Regular analytic classical affine vacuum lift

Read the [formulation](FORMULATION.md), [connected domain](notes/domain.md),
[clock-tube norms](notes/norms.md), [principal lift](notes/lift.md),
[regular lower dictionary](notes/lower.md),
[general-gradient quotient](notes/covariance.md),
[vacuum and limit](notes/vacuum.md) and [scope](notes/scope.md).

With all P8 src directories on sys.path, the independent
native read-only replay is

    python -m p8_affine_vacuum_domain.verify --check

The complete regression is separately run with

    .venv/bin/python -u scripts/p8_snapshot_regression.py

Scientific code, report generation, ordinary replay and CLI
do not monkeypatch SymPy. Only the separately audited full
runner uses its exact GCD adapter. The [report](certificates/regular-analytic-affine-vacuum-lift.json)
pins every own source, proof and test byte.

The regular auxiliary connection has no propagating heavy
spectrum. Neither original P8 nor full V/G/B is closed.
