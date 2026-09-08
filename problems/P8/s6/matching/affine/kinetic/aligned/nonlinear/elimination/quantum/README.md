# First vector-only local quantum coefficient

[S6.47](FORMULATION.md) keeps the temporal constraint and all three
physical vector polarizations in the constant-coefficient determinant.
The [proof](notes/proof.md) derives the pole and finite modified
minimal-subtraction potential, checks the actual clock dependence,
and bounds this local term on the original coefficient tube.

Source alignment does not make the loop determinant a spectator.
An explicit mass-scale choice makes this particular local term small
relative to M^2/tau^2. It does not provide a full curved-background or
quantum remainder, renormalized solution/spectrum, interacting cutoff,
V/G/B conditions, or original P8 closure. No frozen action changes.

Read-only replay with local P8 source roots on PYTHONPATH:

    python -m p8_aligned_quantum.verify --check
