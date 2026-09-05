# S6.2 — First matching routes obstructed

The conventional massive-scalar vacuum realization does **not** yield a
controlled matching of the selected D-only bounce in the tested parent
classes. The [contract](FORMULATION.md) separates these ansatzes from the
broader S6 UV target; [the proof](notes/parent-tests.md) gives their limits.

1. Stable algebraic elimination of affine-X heavy scalars has P_XX>=0.
   The pinned tube has F_XX<0 at every finite time, with F_XX=-10 at the
   bounce in local model units. Matching needs a second-jet correction at
   least as large as the entire target jet, not a small omitted term.
2. A classical minimal Einstein/positive-metric scalar parent obeys
   H_dot<=0 exactly. The bounce requires null stress -4M^2/tau^2; it cannot
   be that parent's physical-metric background solution.
3. A nonminimal scalar-curvature parent with healthy Einstein-frame fields
   can evade the previous obstruction, but not while matching the D Planck
   coefficient with small physical-time C2 errors. The exact target needs
   a common normalized error at least 4/7 across its value, first derivative
   and second derivative. This is necessary, not a sufficient construction.
4. Imposing both nonzero tensor tails excludes even large regular conformal
   or disformal pullbacks of conventional Einstein/NEC scalar parents. The
   auxiliary scale factor cannot grow in both tails while log(a_E) is
   concave. A quantitative three-slice tensor-value test gives necessary
   errors at least 3/5 over t=(-tau,0,tau), or 12/13 over (-2tau,0,2tau).
5. The broader exact map to **quartic Horndeski** also fails globally:
   its necessary Jacobian is C*[1-2/(1+u^2)^3] on the clock, and vanishes at
   u=+-sqrt(2^(1/3)-1). This is failure of the auxiliary map, not a
   pathology in the original bounce. S5's spatial-only map stays regular.

The [frame proof](notes/frame-tests.md) derives both statements and their
different hypotheses. It does not identify all DHOST theories with
Einstein/conventional-scalar parents or exclude arbitrary tube corrections.

Negative controls distinguish a heavy frequency from a stationary Hessian
and exhibit an explicit nonminimal background bounce outside the error
budget. That control has Q->0 tails and is **not** an accepted P8 witness.
No ghost/positivity inference is made from F_XX alone in the actual DHOST
theory. A general UV matching construction remains open.

The [certificate](certificates/parent-tests.json) contains direct lapse-metric
curvature and parent variation, an independent conformal-frame check, exact
heavy-field elimination, the FLINT target sign, tensor-map and necessary
Horndeski-Jacobian derivations, error bounds, counterchecks,
and source/prior pins. Generic analytic lemmas are written, not formalized.
[Primary-source comparison](notes/literature.md) records the known methods.

```sh
PYTHONPATH=problems/P8/src:problems/P8/s5/src:problems/P8/s5/physical/src:problems/P8/s5/control/src:problems/P8/s5/scattering/src:problems/P8/s6/src:problems/P8/s6/matching/src .venv/bin/python -m p8_matching.verify --check
.venv/bin/python -m pytest problems/P8/s6/matching/tests -q
.venv/bin/ruff check problems/P8/s6/matching
git diff --check
```

The verifier only reads; without `--check` it prints a candidate report.
All earlier checkpoints are unchanged. Next research must test a genuinely
different matching mechanism, e.g. leading derivative/curvature interactions
with an explicit spectrum and operator map. Merely adding more conventional
heavy scalars, increasing M*tau or making one of the excluded regular frame
changes does not remove the recorded obstructions.
This is not full S6.2 closure over all parent theories, a D-row exclusion,
a finite-gravity positivity verdict, M1 control or full P8 completion.
