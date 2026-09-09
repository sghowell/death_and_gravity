# Nonminimal conformal scalar

[Formulation](FORMULATION.md) gives the quantitative state-dependent
field inequality, exact cosmological gates, and two different controls.
[Proof](notes/cosmology.md) derives the macroscopic margin with the
field-strength cost retained. [Thermal construction](notes/thermal.md)
supplies an actual short-past scalar SEE realization; it does not
establish the separate future cap.

Run from the repository root after installing the existing P8 environment:

    .venv/bin/python -m p8a_nonminimal.verify --check

The package and its pinned ancestor packages must be on PYTHONPATH;
the repository snapshot runner supplies that setup for regression.
The CLI is read-only; without --check it emits a freshly rebuilt report.
It never writes a certificate or changes scientific SymPy.

The source manifest covers this directory's own Python files, tests,
root Markdown and notes, not child directories or private prototypes.
Original P8 is not closed by this scoped scalar extension.
