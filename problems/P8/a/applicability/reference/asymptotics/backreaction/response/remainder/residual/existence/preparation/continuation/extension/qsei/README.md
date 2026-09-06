# A.15: QSEI on the longer actual semiclassical solution

The same-source A.14 smooth solution now has a freshly proved absolute
QSEI for **all Hadamard target states and proper H2_0 samplers** on its
source-free domain `(x0+L0/2,x0+L)`, with unchanged L0=10^-10 and L=10^-6:

    integral f^2 E_target >= -2 hbar/(16 pi^2) integral |f''|^2.

The actual reference is bounded below by `-hbar/(40pi^2 T0^4)`, rather
than assumed positive there. Weighted metric gains and the new sampler
envelope absorb its possible negative part into the same rounded
coefficient2. The proper-clock, scattering and reference constants are
all independently recalculated; the old short-interval QSEI is not
silently extrapolated.

The available comoving index test still fails on actual curvature. The
sufficient `Q2/duration^2` ratio improves to at least2/5, but that alone
does not supply initial contraction, adequate duration or incompleteness.
P8(a), realistic-field extensions and original P8 remain open.

- [Exact formulation](FORMULATION.md)
- [Signed-reference and two-frequency proof](notes/proof.md)
- [Pinned sources and analytic scope](notes/sources.md)
- [Read-only certificate](certificates/extended-qsei.json)

From the repository root:

```sh
P8_A15=problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/continuation/extension/qsei
.venv/bin/python -m pytest "$P8_A15" -q
.venv/bin/python -m ruff check "$P8_A15"
.venv/bin/python -c 'import runpy,sys; from pathlib import Path; sys.path[:0]=[str(p) for p in Path("problems/P8").rglob("src")]; sys.argv=["p8a_extended_qsei.verify","--check"]; runpy.run_module("p8a_extended_qsei.verify",run_name="__main__")'
```

The report pins A.14 and A.12, replays both full lineages read-only,
and hashes this child's proof, implementation and tests. It is not a
claim that every written functional-analytic step is proof-assistant
formalized or that these field/model constants apply to realistic matter.
