# Source-aligned trace checkpoint

[S6.42](FORMULATION.md) keeps the source-centered mass of S6.41 but
curls W=T-B(phi,x)dphi. The full stationary shifted source has zero
background and first variation, not zero nonlinear value.

The actual quadratic constrained system is exactly the old CD/M1
system plus a positive Proca block. All scalar/vector principal cones
are the original matter cone. The canonical vector equations have a
uniform frequency-squared floor q+1985 for 0<zeta<=1/2000. Explicit
nonlinear source/curl bounds hold on the original closed tube.

See the [source proof](notes/source.md), [constraint proof](notes/dynamics.md)
and [canonical-mode proof](notes/modes.md). These results do not prove
nonlinear secondary constraints, controlled heavy elimination, loop
errors, a cutoff, V/G/B admissibility or original P8 closure.

Read-only replay, with all local source roots on PYTHONPATH:

    python -m p8_affine_aligned.verify --check

The report pins all local sources and fully rebuilds the frozen S6.41
ancestry. Tests include original-variable joint constraints, exact
canonical Euler equations, nonlinear and wrong-sign controls, units,
strict domains and report-mutation rejection. No frozen action changes.
