# Reciprocal-scale geometric obstruction

This new child asks whether the second-metric incompleteness observed in
G1 is specific to its reconstructed lapse. The proof says no
within the explicit reciprocal-scale, positive-shift and source hypotheses
in [FORMULATION.md](FORMULATION.md).

The [proof](notes/proof.md) derives an invariant lapse-order inequality and
an explicit finite affine-length bound. It is not a no-go theorem for the
physical-g bounce or the full P8 operator class. The read-only certificate
is [reciprocal-affine-obstruction.json](certificates/reciprocal-affine-obstruction.json).

Run the tests with:

```
.venv/bin/python -m pytest problems/P8/s6/matching/variable/global/reciprocal/tests -q
```

With the P8 source directories on PYTHONPATH, run
`python -m p8_reciprocal_geometry.verify --check` for the recursive
source-hashed replay. The continuous argument is written and independently
audited, not proof-assistant formalized.
