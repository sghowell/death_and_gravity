# P8 S5 — Perturbative-control extension

S5.1 is implemented: the D-only M0 witness has an exact nonlinear chart
with no lapse derivatives, an independently checked canonical Hamiltonian,
and a lapse-reduced invariant expansion through fourth order. Physical
scale restoration is checked against the covariant background equations.
Local background units compactify both tails and give explicit uniform
bounds for all 34 recorded invariant coefficients, with a uniformly
nonzero normalized lapse Hessian.
The linear checkpoint is unchanged. **S5 as a whole remains open; there is
no strong-coupling or UV verdict yet.**

Read [the acceptance contract](FORMULATION.md),
[the derivation and its limitations](notes/nonlinear-d.md),
[the audit record](notes/audit.md), and
[the primary-source method comparison](notes/literature.md).

## Reproduce

From the repository root, using the existing locked environment:

```sh
.venv/bin/python -m pytest problems/P8/s5/tests -q
PYTHONPATH=problems/P8/src:problems/P8/s5/src .venv/bin/python -m p8_s5.verify --check
PYTHONPATH=problems/P8/src .venv/bin/python -m p8.verify_all --check
.venv/bin/ruff check problems/P8/s5
git diff --check
```

The [S5.1 report](certificates/d-nonlinear.json) stores exact coefficients,
negative-control-tested identity checks, all-time denominator/sign proofs,
source hashes and the old classification's hash. The verifier only reads:
running it without `--check` prints a proposed report for inspection, not an
automatic overwrite. Existing P8 tests and certificate replays remain usable.

## What comes next

The invariant generator retains full spatial curvature, the volume factor,
and trace/traceless momentum components. To obtain physical interaction
vertices, expand that geometry, solve the spatial momentum constraints,
and construct canonical scalar/tensor modes in regular overlapping charts.
Then evaluate cubic exchange plus quartic contact contributions and prove
appropriate tail bounds. Do not estimate the cutoff from the stored
off-shell coefficients or from the lapse Hessian alone.

The useful exact checks are h_NN=-2J with J>0 at all finite times,
h_NN(0)=-4, and 1/(2d)<J<=2/d. The unweighted Hessian vanishes in the tails,
but local background units give -4<=h_local,NN<-1 and a uniform local
algebraic lapse branch. This is not a physical-mode cutoff. A nonzero background
reference is E_ref=1/(tau sqrt(d)); the curvature scale stays between
sqrt(2) and sqrt(6) times E_ref. M tau is free, but no quantitative cutoff
hierarchy has been established by that freedom alone.
