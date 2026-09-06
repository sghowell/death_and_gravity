# P8(a) A.11 — smooth local SEE through conserved joint preparation

The full actual semiclassical map has a certified contraction on a
dimensionless conformal slab of length `10^-10`, starting at the old
prepared plateau label `y=5/2`. A smooth conserved preparation source
turns off halfway through. The final open half-slab solves the **full
unforced SEE**, with the original vacuum jointly transported on the new
smooth metric and the ordinary radiation normalization unchanged.

The main quantitative caps are full-map contraction `<3/25`, center
image `<63*10^-9`, and a self-map inside the radius-`10^-6` potential-
derivative ball. A common flat-start Picard argument proves smoothness
on the same slab; the generic logarithmic inverse is not asserted to
preserve `C1`. The positive inverse pole is retained.

- [Precise theorem and physical scope](FORMULATION.md)
- [Actual map, conserved source and energy constraint](notes/construction.md)
- [Complete rational contraction bounds](notes/contraction.md)
- [All-order smoothness argument](notes/regularity.md)
- [Read-only replay certificate](certificates/smooth-preparation.json)
- [Primary-source audit](notes/sources.md)
- [Remaining P8(a) obligations](notes/remaining.md)

From the repository root:

```sh
.venv/bin/python -m pytest problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/tests -q
.venv/bin/python -c 'import runpy,sys;from pathlib import Path;sys.path[:0]=[str(p) for p in Path("problems/P8").rglob("src") if p.is_dir()];sys.argv=["p8a_preparation.verify","--check"];runpy.run_module("p8a_preparation.verify",run_name="__main__")'
```

Without `--check`, the verifier prints the reproducible report; it never
overwrites the checked-in evidence. Source hashes cover the implementation,
tests and written proof, and the immutable A.10/A.8/A.7 lineage is replayed.

This is not whole-window continuation, automatic transfer of A.9's QSEI,
a realistic-field focusing theorem, or closure of original P8. The
preparation stress need not obey energy conditions and is not globally
compact in the unchanged off-shell past.
