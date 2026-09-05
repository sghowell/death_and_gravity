# P8 S5.3b / S5.4 — Hard finite-time tree control and UV applicability

The D-only M0 witness meets a **specified finite-time, fixed-total-momentum,
nonexceptional hard-channel tree criterion through quartic order**, uniformly
at every background time. The proof includes scalar/tensor cubic exchange
and quartic Hamiltonian contacts, with exact free mode evolution.

An intentionally loose sufficient choice is **M tau >= 10^270**. At each
time, the band `10^8<=ell0*k_physical<=10^9`, with the stated nonzero-transfer
restriction, admits windows of half-length `ell0/100` whose specified
connected cubic/quartic tree transition blocks have norm <=`10^-3`.
This is an existence-grade sufficient condition, not a realistic parameter
estimate, a necessary scale, inclusive scattering or an all-orders cutoff.

Read [the operational contract](FORMULATION.md),
[the finite-time proof](notes/finite-time.md),
[the separate frozen-symbol calculation](notes/frequency.md), and
[the UV applicability audit](notes/uv-applicability.md).
The [report](certificates/tree-control.json) stores exact majorants, five
frozen tree regressions, a free-field cancellation control and source pins.

The S5.4 applicability audit is complete, but positivity is **not applied**:
the smooth clock-tube contract supplies neither a stationary vacuum/scattering
matching problem nor gravitational dispersion assumptions. A UV classification
requires an additional acceptance contract. Inclusive/forward and homogeneous
channels, vacuum production/full Fock norms, global evolution/vacuum choice,
  loops/all orders, M1 interaction control, nonlinear stability and the full
P8 program remain open. The older pinned checkpoints are unchanged.

## Reproduce

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src .venv/bin/python -m p8_scattering.verify --check
.venv/bin/python -m pytest problems/P8/s5/scattering/tests -q
.venv/bin/ruff check problems/P8
git diff --check
```

The verifier only reads; without `--check` it prints a candidate report.
Direct frequency seeds are checked against multilinear expansion of the
independently implemented S5.2 velocity kernels.
