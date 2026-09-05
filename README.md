# death_and_gravity

AI-assisted attack on ten open problems in theoretical cosmology, organised
around a pre-registered problems document and a certificate-gated ledger.

- Problems document: `docs/problems/open-problems-theoretical-cosmology-2026.tex`
  (frozen at tag `problems-v1.1`; every later change goes in its Revision log).
- Starting-point assessment and pre-work: `docs/assessment-2026-08-21-starting-points.md`.
- Claims ledger: `CLAIMS.md` — nothing is stated as a result unless it is
  recorded there with the evidence path and a level.
- P8(a): [conditional local-QEI focusing optimization](problems/P8/a/README.md).
  Exact localization and a certified two-parameter trial lower the same-partition
  cubic contraction threshold by 8.825--8.826%. Local QEI validity duration,
  initial Ricci assumptions and physical field/state calibration remain explicit;
  this is not full P8(a) or a cosmological application.
  The [exact radiation-background difference QSEI](problems/P8/a/applicability/README.md)
  proves a Minkowski-coefficient bound on every compact positive-time sampling
  interval in that fixed testbed, without assuming a short averaging duration.
  The new [reference-stress and absolute-QSEI calculation](problems/P8/a/applicability/reference/README.md)
  computes the positive reference effective energy and removes the conditional
  reference-stress parameter for this testbed. Its reference state is not an
  exact semiclassical source for the radiation metric; self-consistent
  cosmological application remains open.
- P8(b): [completed scoped classification and verification commands](problems/P8/README.md).
  All 32 M0/M1 sector-row verdicts are decided: minimal optional groups are
  C or D without extra matter, and CD with rolling canonical matter.
  This is an all-time linear-principal result, not UV or nonlinear completion.
- P8 S5: [nonlinear D-witness preparation](problems/P8/s5/README.md),
  [physical cubic/quartic scalar-tensor reduction](problems/P8/s5/physical/README.md),
  [exact oscillator/adiabatic-indicator bounds](problems/P8/s5/control/README.md),
  and [hard finite-time tree control / UV applicability](problems/P8/s5/scattering/README.md).
  The D-only finite-time, nonexceptional cubic/quartic tree criterion now has
  a uniform sufficient scale choice. Inclusive scattering, all-orders control,
  M1 interaction control and UV admissibility remain open; no UV completion verdict.
- P8 CD/M1: [nonlinear preparation and principal-cone robustness](problems/P8/s5/matter/README.md).
  The matter branch now has a regular nonlinear spatial chart, its required
  canonical boundary, seven-velocity Legendre reduction and quartic invariant
  coefficients with all-time bounds. Its new
  [physical cubic/quartic reduction](problems/P8/s5/matter/physical/README.md)
  solves all three matter-sourced spatial constraints and retains the paired
  scalar-matter canonical shifts. The unnormalized velocity kernels include the
  full coupled scalar Legendre contact and both tensor polarizations. The new
  [coupled canonical normalization and free-energy bound](problems/P8/s5/matter/physical/control/README.md)
  retains finite-momentum mixing and all time-dependent generators, and gives
  a uniform sufficient local band with energy variation below a factor of two.
  Its conservative q>=10^20 threshold is not an interacting cutoff;
  M1 interaction control remains open.
  Exact mixing mismatch forces a superluminal principal mode with frozen
  canonical matter; a specified one-sided covariant deformation stays causal
  within the exceptional relation. Neither result computes loop corrections.
- P8 S6: [conditional UV framework and vacuum-extension test](problems/P8/s6/README.md).
  A smooth clock tube does not fix the vacuum positivity coefficient. The
  adopted vacuum/Regge framework requires controlled matching to the bounce;
  a healthy spliced vacuum is not a UV-completion result.
  [First parent-matching tests](problems/P8/s6/matching/README.md) obstruct
  conventional scalar parents, including regular frame changes with the
  required tensor tails. An exact D-to-quartic-Horndeski map loses
  invertibility at two finite times; this is not a singularity of D itself.
  General derivative/curvature or quantum matching mechanisms remain open.

## Layout

```
docs/                       problems document, assessment, revision history
problems/P<n>/              one directory per problem actually being worked
  FORMULATION.md            frozen statement of the entry point being attacked
  notes/                    literature digests, derivations
  data/                     frozen inputs (+ MANIFEST.json with sha256, URLs, upstream commits)
  src/                      code
  certificates/             machine-checkable artifacts backing ledger entries
CLAIMS.md                   the ledger
```

## Ledger levels

Borrowed from the Empiricist harness (`~/dev/empiricist`); promotion requires
machine evidence at every step above HEURISTIC.

| Level | Meaning |
|---|---|
| HEURISTIC | numerical/analytic lore, no certificate |
| CONJECTURED | precise statement, tested numerically, not certified |
| VERIFIED_N | independently re-run by a second implementation/engine |
| CERTIFIED | machine-checkable certificate (interval enclosure, exact dual, exhaustive enumeration) committed under `certificates/` |
| FORMALIZED | kernel-checked proof (Lean 4) |
| REFUTED | terminal |

## Building the problems document

```
brew install tectonic
cd docs/problems && tectonic open-problems-theoretical-cosmology-2026.tex
```

## Python environment (P9)

```
uv sync
uv run python -c "import cvxpy, flint; print(cvxpy.__version__)"
```
