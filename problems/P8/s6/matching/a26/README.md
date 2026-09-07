# S6.19.A26 — literal constant-clock domain audit

The June-2026 two-field DHOST action has a generic nonquadratic
constant-clock directional limit, not merely individual coefficient poles.
Full symmetric-Hessian variation and a separate Fourier/IBP bulk check show
that no boundary term restores a bilinear vacuum Hessian. The literal
benchmark cannot remain unchanged at nonzero X near zero while acquiring a smooth
vacuum chart on an open finite-phi neighborhood.

This does not exclude a separately named smooth off-tube repair or certify
a cosmological instability, positivity failure, or UV no-completion. The
fractional Q/W2 note retains the complete reconstructed-source cancellations.

Read [FORMULATION.md](FORMULATION.md), [proof](notes/proof.md), and
[primary-source audit](notes/sources.md). The report is
[a26-vacuum-domain.json](certificates/a26-vacuum-domain.json).

From the repository root:

```sh
.venv/bin/pytest -q problems/P8/s6/matching/a26/tests
.venv/bin/python -c 'import pathlib,runpy,sys; sys.path[:0]=[str(p) for p in pathlib.Path("problems/P8").rglob("src")]; runpy.run_module("p8_a26_vacuum.verify",run_name="__main__")' --check
.venv/bin/ruff check problems/P8/s6/matching/a26
```

The verifier only reads inputs and writes its result to stdout. Without
`--check` it prints a candidate report; it never regenerates a certificate
on disk. S6.1 and the adopted S6 contract are pinned/replayed, and all direct
new sources, tests and proof notes are hashed. New descendants are not
silently included in this checkpoint's source scope.
