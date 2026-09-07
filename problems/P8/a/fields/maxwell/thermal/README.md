# A.19 — actual thermal-photon semiclassical history

An actual positive Hadamard state of the two physical photon
polarizations realizes an exact semiclassical history inside A.18's
anchored radiation C3 tube. This optional example uses the explicitly
named beta_M=0, Lambda=0 prescription, with no added classical matter.
It does not change the independent all-Hadamard incompleteness theorem
or reopen that scoped objective.

[The formulation](FORMULATION.md) states the exact model and quantitative
history bounds. [The proof](notes/proof.md) constructs the physical state,
fixes its thermal charge from the Einstein constraint, derives both
independent SEE components and controls the entire past interval. It
also identifies the actual smooth open low branch and its finite proper
endpoint. The physical anomaly is retained; this is not an order-reduced
or formal approximate solution.

The thermodynamic b_T is a conformal-time length, not the photon finite
coefficient beta_M or an older scalar gamma. Its physical temperature
obeys k_B T=hbar/(a b_T), with c=1. The physical density is Q/a^4 plus
the nonzero vacuum anomaly. [The source audit](notes/sources.md) records
the field-strength construction and source/sign/temperature dictionary.

## Read-only verification

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/a/fields/maxwell/thermal/tests
.venv/bin/ruff check problems/P8/a/fields/maxwell/thermal
.venv/bin/python -c 'from pathlib import Path; import sys; root=Path("problems/P8/a/fields/maxwell/thermal"); m=root.parent; f=m/"focusing"; c=f/"cosmology"; old=m.parents[1]/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"; sys.path[:0]=[*(str(p/"src") for p in (root,m,f,c)),str(old.parents[2]/"qsei"/"src"),*(str(p/"src") for p in (old,*old.parents))]; from p8a_maxwell_thermal.verify import main; main()' --check
```

The CLI compares [the frozen report](certificates/thermal-see.json) with
a fresh exact replay, validates the A.18 report and full A.17/A.16
lineage, and checks this child's source hashes without writing files.
Without `--check`, it prints a candidate only: it never repairs or
overwrites the report. The independent Fraction reconstruction and
parent-owned covariant audit are part of the pinned source inventory.

## Interfaces and limits

- `state.prescription(beta_m=0, cosmological_constant=0)` explicitly
  admits only this example's finite prescription and source choice.
- `state.thermal_amplitude(hbar=..., b_T=...)` computes the physical
  two-polarization thermal charge for exact positive rational inputs.
  The written state construction permits any positive real b_T, including
  the exact fourth root fixed by the constraint; this numerical helper's
  rational domain does not restrict that theorem.
- `dynamics.branch_point(y, delta)` supplies exact interior metric/stress
  data and rejects the critical or high algebraic branch.
- `dynamics.clock` and `jet_functions` are symbolic exact identities;
  `endpoint_data(delta)` gives the exact primitive and analytic rational
  endpoint interval, not a floating-point endpoint scan.
- `bounds.history(delta)` gives continuous C3 enclosures for every
  x in [-1/100,0] and every 0<delta<=10^-8.

Numerical interfaces reject bools, binary floats, nonfinite values and
unresolved symbols. Formal stress/clock/identity expressions remain
symbolic. The finite endpoint is excluded from the regular field/metric
domain; no sampling through it or smooth branch switch is claimed.
This particular solution has positive EED, so its explicit endpoint is
not a new QEI-only proof. Fundamental EFT validity at the high-curvature
endpoint, generic-beta actual existence, interacting QED and observed
cosmological fitting remain outside this optional example.
