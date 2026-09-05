# P8 S5.3 — Control test, underway

S5.3a now derives the exact scalar/tensor oscillator frequencies after
including both kinetic normalization and the expanding `a^3` volume factor.
It proves conservative, uniform bounds on the first two local adiabatic
indicators on covering scalar charts and for tensors, including both compact
tail limits. **This does not establish interaction control or a cutoff.**

The generic [physical cubic/quartic kernels](../physical/README.md) are
available. The next required work is their local frequency-space assembly,
cubic-exchange plus quartic-contact amplitudes, and a frequency/angular/IR
domain with justified interaction and adiabatic errors.

Read [the oscillator derivation and bound proof](notes/oscillators.md).
The [exact report](certificates/adiabatic.json) pins the earlier physical
interaction report and all of this stage's code, tests and notes.

## Reproduce

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src .venv/bin/python -m p8_control.verify --check
.venv/bin/python -m pytest problems/P8/s5/control/tests -q
.venv/bin/ruff check problems/P8/s5/control
```

The present sufficient threshold `q>=10^14`, equivalently
`k_physical>=10^7/ell`, is an intentionally loose polynomial-majorant bound,
**not** a necessary adiabatic threshold or an EFT cutoff. It only establishes
small linear-mode indicators at those local frequencies. It says nothing
yet about whether those frequencies lie below the interacting theory's
cutoff; proving a nonempty overlap is precisely the unfinished S5.3 task.
