# death_and_gravity

AI-assisted attack on ten open problems in theoretical cosmology, organised
around a pre-registered problems document and a certificate-gated ledger.

- Problems document: `docs/problems/open-problems-theoretical-cosmology-2026.tex`
  (frozen at tag `problems-v1.1`; every later change goes in its Revision log).
- Starting-point assessment and pre-work: `docs/assessment-2026-08-21-starting-points.md`.
- Claims ledger: `CLAIMS.md` — nothing is stated as a result unless it is
  recorded there with the evidence path and a level.
- P8(b): [completed scoped classification and verification commands](problems/P8/README.md).
  All 32 M0/M1 sector-row verdicts are decided: minimal optional groups are
  C or D without extra matter, and CD with rolling canonical matter.
  This is an all-time linear-principal result, not UV or nonlinear completion.

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
