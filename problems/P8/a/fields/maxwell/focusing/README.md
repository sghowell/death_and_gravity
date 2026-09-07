# A.17 — photon QSEI with a geometric contraction history

A field-specific conditional incompleteness theorem with no initial
pointwise SEC/Ricci-sign premise. Its actual geometric hypotheses are a
quantified short contraction history and relative future Hubble-jet bounds
on any normal segment which reaches the proposed endpoint. The bounds allow
inverse-power growth near that endpoint, and the exact quantum constants
are compatible with macroscopic sampling scales.

[The formulation](FORMULATION.md) states every quantifier, the real beta_M
prescription family, actual SEE and signed extra-source requirements.
[The proof](notes/proof.md) gives the two-sided sampler identity, global
FLRW contradiction and a smooth future-complete geometric countercontrol.
[The source audit](notes/sources.md) separates the new result from the
initial pointwise hypotheses and compact-Cauchy theorems in earlier work.

The named rational calibration has zero-source margin >3/4 and still >1/3
with the displayed extra-source allowance. It is not an observed-universe
or thermal-radiation witness. The complete geometric countercontrol is not
a permitted Maxwell SEE solution. No larger-spacetime inextendibility,
fundamental EFT validity or original P8 closure is claimed.

## Read-only verification

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/a/fields/maxwell/focusing/tests
.venv/bin/ruff check problems/P8/a/fields/maxwell/focusing
.venv/bin/python -c 'from pathlib import Path; import sys; root=Path("problems/P8/a/fields/maxwell/focusing"); m=root.parent; old=m.parents[1]/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"; sys.path[:0]=[str(root/"src"),str(m/"src"),str(old.parents[2]/"qsei"/"src"),*(str(p/"src") for p in (old,*old.parents))]; from p8a_maxwell_focusing.verify import main; main()' --check
```

The CLI checks [the frozen report](certificates/history-focusing.json), the
A.16 report hash and its full source-hashed lineage. It replays exact
symbolic identities, independent Fraction integration and all new source
hashes without writing anything. Without `--check` it prints a candidate
report, not a file repair. The separate parent-owned covariant audit is
included unchanged.

## Calculation interfaces

- `history.moments`, `history.lower`: exact cubic moments and the quantified
  history lower bound, not an inference from instantaneous K alone.
- `envelope.coefficients`: separate past/future jet caps and explicit beta_M;
  returns the relative-pole and reference-anomaly costs.
- `focusing.theorem_constants`: a rational sufficient-test margin with delta
  and the separately supplied source sigma.
- `focusing.physical_enclosure`: the physical signed-source dictionary and
  rational pi>3 upper bound on delta; no hidden choice of units or sources.
- `controls.calibration`: exact smoothing/jet and source-mismatch constants
  for the future-complete geometric countercontrol.

Exact numerical interfaces reject booleans, binary floats, nonfinite or
unresolved symbolic quantities. A positive returned margin does not verify
that a supplied metric, state or source satisfies the theorem hypotheses.
