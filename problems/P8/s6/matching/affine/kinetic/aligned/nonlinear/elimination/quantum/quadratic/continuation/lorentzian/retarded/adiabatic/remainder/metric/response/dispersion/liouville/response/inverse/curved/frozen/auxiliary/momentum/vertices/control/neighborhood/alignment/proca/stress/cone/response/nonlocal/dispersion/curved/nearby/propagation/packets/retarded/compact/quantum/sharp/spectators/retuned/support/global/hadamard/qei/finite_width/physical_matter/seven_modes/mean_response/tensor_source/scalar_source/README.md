# Actual coupled scalar mean-source certificate

This increment derives the remaining two-mode scalar sources
from the actual action before homogeneous mean restriction.
Both canonical shifts, the time boundary, all four source
terms and the physical matter fluctuation terms are retained.

Read the [formulation](FORMULATION.md), [action](notes/action.md),
[phase and state](notes/phase.md), [mean](notes/mean.md),
[observable](notes/observable.md), [bounds](notes/bounds.md),
[center](notes/center.md) and [scope](notes/scope.md).

With the repository environment and all P8 src directories
on sys.path, the read-only native replay is

    python -m p8_scalar_mean_source.verify --check

The separate complete regression command is

    .venv/bin/python -u scripts/p8_snapshot_regression.py

Only that separately audited full-regression runner uses its
exact GCD adapter. Native scientific checks, report generation,
ordinary replay and CLI do not monkeypatch SymPy.
Every own source, proof, formulation and test byte is pinned.

The [report](certificates/actual-coupled-scalar-mean-response.json)
contains exact source-kernel fingerprints, not a replacement
of the kernels by a finite-momentum approximation. Kernels
are fully rebuilt before hashing. This is controlled leading
response near the bounce, not a uniform-tail or UV theorem.
Original P8 remains open.
