# P8 S5.2 — Physical scalar/tensor interactions

The D-only M0 witness now has an exact, reproducible generator of reduced
cubic and quartic Hamiltonian and velocity kernels. It solves all three
spatial momentum constraints, keeps generic non-collinear spatial momenta
and both tensor polarizations, and includes the time-dependent canonical
boundaries, normalization derivatives and quartic Legendre correction.

This completes the **nonexceptional-mode interaction reduction**, not S5's
perturbative-control or UV gates. Zero-momentum channels, cubic exchange,
on-shell amplitudes, adiabatic control and a quantitative all-time cutoff
hierarchy remain open. Fixed-momentum kernel bounds are not such a hierarchy.

Read [the precise scope](FORMULATION.md), [the derivation](notes/reduction.md)
and [the audit/next steps](notes/audit.md). The
[report](certificates/interactions.json) records 14 exact symbolic-time phase
kernels with conservative all-time bounds and 10 canonically normalized
velocity examples, with separate scalar/tensor quartic Legendre terms.
These examples are regressions of a generic procedural construction, not
an exhaustive sampling argument for all interactions or scattering angles.

## Reproduce

From the repository root:

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src .venv/bin/python -m p8_physical.verify --check
.venv/bin/python -m pytest problems/P8/s5/physical/tests -q
PYTHONPATH=problems/P8/src:problems/P8/s5/src .venv/bin/python -m p8_s5.verify --check
PYTHONPATH=problems/P8/src .venv/bin/python -m p8.verify_all --check
.venv/bin/ruff check problems/P8
git diff --check
```

Running the verifier without `--check` prints a candidate JSON report;
it does not write files. The report pins this source tree, the published
S5.1 report and its inputs, and the earlier classification's inputs.
This nested directory deliberately preserves both published source-glob
checkpoints. Their historical descriptions of work then remaining are
unchanged; this directory records the subsequent reduction.

## Example API

With the same `PYTHONPATH`:

```python
from p8_physical.vertices import Leg
from p8_physical.lagrangian import kernel

legs = tuple(Leg(k, 's_dot') for k in
             ((10, 20, 0), (0, 10, 30), (-10, -30, -30)))
result = kernel(legs, time_point=0, chart='gamma', normalize=True)
print(result['kernel'])
```

Coordinates are `s`, `t`; velocities `s_dot`, `t_dot`. Phase kernels also
accept scalar `p` and tensor `pi` momenta. Supply a rational TT matrix for
each tensor leg; `tensor_basis(wave)` supplies two possible matrices.
`time_point=None` retains exact compact time `x`; otherwise supply an exact
rational in `[-1,1]`. A returned coefficient is multilinear in **labelled**
external legs, with no extra `1/n!`. This convention retains repeated-field
combinatorial factors and must be respected when building amplitudes.
