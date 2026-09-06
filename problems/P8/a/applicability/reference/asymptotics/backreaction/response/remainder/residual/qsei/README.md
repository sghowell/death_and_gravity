# P8(a) A.9 — All-sampler QSEI on the actual prepared metric

For the unchanged A.8 metric and named prescription at delta<=10^-14,
every Hadamard target state and every real compact sampler in 2<y<3 obey

    integral h² E_omega dt >= -5*hbar/(16*pi²)*integral |h_ddot|² dt.

The proof controls the actual two-point spectral functional in both the
infrared and ultraviolet. It uses only two sampler derivatives; the
target state need not share the reference state's symmetries. The
constant is sufficient, not sharp. No exact SEE solution or new
cosmological incompleteness theorem is claimed.

Read the [fixed contract](FORMULATION.md), [proof](notes/proof.md) and
[source boundary](notes/sources.md). The [certificate](certificates/prepared-qsei.json)
pins all local source and audit files plus A.8's unchanged certificate.

## Reproduce

From the repository root:

```sh
.venv/bin/python -m pytest problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/qsei/tests -q
.venv/bin/ruff check problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/qsei
```

For the read-only CLI, expose the nine A-track source directories:

```sh
PYTHONPATH=problems/P8/a/src:problems/P8/a/applicability/src:problems/P8/a/applicability/reference/src:problems/P8/a/applicability/reference/asymptotics/src:problems/P8/a/applicability/reference/asymptotics/backreaction/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/src:problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/qsei/src .venv/bin/python -m p8a_qsei.verify --check
```

Without `--check` the verifier prints a candidate report and never writes
one. The underlying Volterra and positive-type analytic lemmas are stated
in the proof rather than attributed to finite test samples.
