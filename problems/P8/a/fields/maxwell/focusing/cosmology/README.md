# A.18 — robust photon cosmological-strength calibration

The A.17 incompleteness theorem now has an explicit C3 neighborhood of
radiation/matter-era geometric histories with a sourced margin greater
than 1/8. Every finite photon beta_M remains explicit through
delta(1+|beta_M|)<=10^-8. This is a free-photon/global-flat-FLRW theorem
specialization with cosmological-scale constants, not an observed-universe
fit or a new exact quantum solution.

[The formulation](FORMULATION.md) states the actual past tube, separate
conditional future envelope, signed cosmological/additional-source budget
and proper-time conclusion. [The proof](notes/proof.md) treats the whole
continuous power interval and gives smooth future-complete comparison
metrics satisfying the same geometric bounds. They are not allowed
quantum SEE witnesses: their role is to prevent a geometric contradiction
from being mislabeled as a photon-QSEI result. The rejected tighter future
envelope is preserved as an exact negative control.

The reference H_* tau=2 is distinguished from the actual unanchored
observer's H_0 tau in [199/100,201/100]. The source and finite prescription
are not silently changed. [The source audit](notes/sources.md) maps the
unchanged field/index inputs, and [the closure assessment](notes/closure.md)
states the precise relation to the literal original objective without
assigning root P8 completion.

## Read-only verification

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/a/fields/maxwell/focusing/cosmology/tests
.venv/bin/ruff check problems/P8/a/fields/maxwell/focusing/cosmology
.venv/bin/python -c 'from pathlib import Path; import sys; root=Path("problems/P8/a/fields/maxwell/focusing/cosmology"); f=root.parent; m=f.parent; old=m.parents[1]/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"; sys.path[:0]=[str(root/"src"),str(f/"src"),str(m/"src"),str(old.parents[2]/"qsei"/"src"),*(str(p/"src") for p in (old,*old.parents))]; from p8a_maxwell_cosmology.verify import main; main()' --check
```

The CLI compares [the frozen report](certificates/cosmological-calibration.json)
with a fresh exact replay, checks the A.17 SHA and its full source-hashed
lineage, and checks this child's source hashes without writing files.
Without `--check` it only prints a candidate report; it does not repair
or regenerate a certificate on disk. The independent Fraction engine and
parent-owned covariant audit are included in the hashed sources.

## Calculation interfaces

- `geometry.history_point`, `geometry.neighborhood`: exact all-p history
  and robust C3 margins, with an optional explicitly anchored observer.
- `calibration.affine_cost(past_caps, future_caps)`: separate past/future
  caps, returning C0 and Cbeta from the immutable photon envelope.
- `calibration.theorem_gate(weighted_delta, sigma)`: rational sufficient
  margin; `weighted_delta` bounds delta(1+|beta_M|), not delta alone.
- `dictionary.proper_scales`, `dictionary.source_budget`: actual versus
  reference clocks, explicit signed source loss and rational pi>3 gate.
- `controls.complete_geometry`, `controls.rejected_tight_future`: complete
  geometric comparisons and the exact geometry-only exclusion control.

Exact numerical inputs reject bools, binary floats, nonfinite values and
unresolved symbols. Formal identity/switch functions remain symbolic.
No positive returned margin verifies an actual supplied metric, state or
future-extension hypothesis by itself. The conditional theorem is the
analytic implication proved in the notes, not a claim of all-field,
all-spacetime, optimal-constant or fundamental EFT validity.
