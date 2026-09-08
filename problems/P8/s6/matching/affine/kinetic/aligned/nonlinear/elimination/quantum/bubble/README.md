# Vector two-insertion quantum component

[Scope](FORMULATION.md), [pole derivation](notes/pole.md),
[finite remainder proof](notes/remainder.md).

S6.48 retains the full Proca numerator in the actual S6.42 mass
direction. It gives the momentum-dependent quadratic pole and a
finite loop-only bound after a specified local Taylor subtraction.
It does not give a full curved-background quantum correction or
close original P8. No frozen ancestor or counterterm is changed.

From the repository root, the ordinary suite is
`PYTHONHASHSEED=0 .venv/bin/python -m pytest
problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/bubble/tests
-q -p no:faulthandler`. The test conftest discovers source roots.
For standalone read-only replay, expose the P8 source roots and run
`python -m p8_vector_bubble.verify --check`. The CLI never writes.
