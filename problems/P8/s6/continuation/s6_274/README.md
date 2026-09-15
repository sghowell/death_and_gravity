# P8 S6.274: the same finite hybrid's core and regulator comparisons

This packet improves estimates for the **same four S6.273 finite hybrid
solutions**. It does not define a new Hamiltonian, state, physical cutoff
or background history.

On the unchanged interval |u| <= 1e-180, the evolved positive coherent
outside-core probability is < 1e-6. At fixed ordering, changing between
the original c=1,2 cutoffs changes the five live homogeneous coordinates,
phase-factored state and physical volume by < 1e-580, 1e-12, 1e-514,
respectively. At fixed cutoff, changing between the original calibrated
coherent and Weyl orderings gives < 1e-630, 1e-65, 1e-566.

The complete constrained quadratic energy is derived from the ORIGINAL
Hamiltonian, including trace/shear cancellation, the full lapse response,
every primitive/heavy/Gauss contact and complete tensor curvature.
A larger COMPLEX ANALYSIS neighborhood controls the whole nonlinear
remainder. It does not enlarge the physical cutoff. An exact sixth
antiheat identity retains its full positive-heat remainder and supplies
stronger bounds for both unchanged operators.

Original V/G/B/P8 remain OPEN. These are positive-POVM and finite-regulator
comparisons, not exact support, regulator removal, homogeneous
quantization, original fully quantum dynamics, matching or global closure.
The S6.273 turnaround is not upgraded to a unique or strict minimum.

## Contents and replay

The19 sources are8 package modules,2 test files, this README, FORMULATION
and7 proof notes. The report has20 top-level fields,61 named exact checks,
289 scalar entries,57 gates,13 controls including496 rejected inputs,
9 unchanged primitive frontier rows and130 distinct matching rows.
The independent science suite has1210 tests; the complete packet adds25
read-only certificate tests.

From the repository, with its existing environment and replay harness:

```sh
.venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_274
.venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_hybrid_core_regulator_comparison.verify
```

Native/direct/ordinary/CLI use original SymPy. The exact-GCD adapter is
restricted to a separately captured FULL regression. Reports and science
sources are immutable after freezing; a correction requires a successor
packet or external assessment, not an in-place rewrite.

Start with FORMULATION.md, then notes/quadratic.md and notes/operators.md
for the two principal arguments. notes/dynamics.md explains the coupled
comparison and the necessary scalar-phase convention.
