# A.8: finite actual semiclassical residual on the prepared radiation metric

This checkpoint bounds the actual renormalized scalar stress against the
positive radiation reference and derives a finite-amplitude residual of
the actual Einstein-plus-radiation equation on the unchanged prepared
metric. It uses the exact conformal clock, conserved trace anomaly and
specified past state; no numerical response input is assumed.

In the declared prescription, `delta<=10^-14` gives less than one-percent
error in each positive reference stress component. The corresponding
residual is below `10^-17` of each specified positive zeroth-order
classical component. Physical epsilon=1 is included when
`t_star>=2*10^7*sqrt(d)`. These are fixed-testbed error bounds, not an
exact solution, stability result or cosmological theorem.

Read [FORMULATION.md](FORMULATION.md), [the proof](notes/proof.md), and
[the source boundary](notes/sources.md). Earlier files remain unchanged.

```sh
PYTHONPATH=problems/P8/a/src:problems/P8/a/applicability/src:problems/P8/a/applicability/reference/src:problems/P8/a/applicability/reference/asymptotics/src:problems/P8/a/applicability/reference/asymptotics/backreaction/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/src .venv/bin/python -m p8a_residual.verify --check
.venv/bin/python -m pytest problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/tests -q
```

The verifier prints a report without writing when run without `--check`.
The pinned report is `certificates/radiation-residual.json`. Generic symbolic
reconstruction retains finite terms and the past integration constant;
the numerical theorem applies only to its explicitly named model and
preparation. QSEI, exact/nearby SEE control and P8(a) completion remain open.
