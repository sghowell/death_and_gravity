# Composite all-mode principal-cone screen

Under the explicit positive-affine, regular common-flat assumptions,
requiring both tensor modes to be no faster than composite matter and
the relative vector to have strictly positive principal speed forces
the physical Hubble to be nonincreasing. This excludes any degeneracy
of contraction-to-expansion transition within that stated class.

See [FORMULATION.md](FORMULATION.md), [proof](notes/proof.md) and
[certificate](certificates/composite-cones.json). The exact source/clock
identities and independent Fraction controls are reproducible with:

```sh
.venv/bin/python -m pytest problems/P8/s6/matching/composite/perturbations/cones -q
.venv/bin/python -c 'import runpy,sys;from pathlib import Path;sys.path[:0]=[str(p) for p in Path("problems/P8").rglob("src") if p.is_dir()];sys.argv=["p8_composite_cones.verify","--check"];runpy.run_module("p8_composite_cones.verify",run_name="__main__")'
```

All certificate commands are read-only. Existing source-hashed ancestors
are immutable. No general composite/UV exclusion, finite-band matching
verdict or completion of original P8 is asserted.
