# S6.5.HR: root-covered normalized-Hubble monotonicity

The same regular common-flat HR parent class as S6.4 cannot pass from
physical contraction to expansion, even through an arbitrarily degenerate
zero or a static plateau. The result holds on any connected regular
interval and sharpens the CD mismatch to an endpoint-H bound alone.
It is a new child certificate; S6.4 and all earlier actions remain unchanged.

## Exact hypotheses

Use the pinned HR action, `+---` physical matter metric `g`, positive
Einstein coefficients `G=M_g^2,F=M_f^2`, positive `m4`, and arbitrary
constant real `beta0,...,beta4`. Both metrics have a common spatially flat
FLRW foliation with smooth, positive, finite scale factors `a,b` and lapses
`N_g,N_f` at every point of the connected interval. The positive square-root
branch has eigenvalues `(c,y,y,y)`, with `c=N_f/N_g>0`, `y=b/a>0`.

The primary matter content is separately conserved, minimally coupled
classical NEC matter on `g` only. The extension permits independent
minimally coupled, separately conserved matter sectors on each metric,
each obeying its own NEC. This is not a shared field coupled to both
metrics. The physical geometry and proper time are not redefined.

## Monotonicity statement

Let `P(y)=m4(beta1+2 beta2 y+beta3 y^2)` and
`D_g=N_g^{-1}d/dt`.

- If `P` is not the zero polynomial, define
  `A=G+F y^2` and `Z=H_g/sqrt(A)`. Then `D_g Z<=0` everywhere on the
  connected interval. On the open dynamical set `P(y)!=0`,
  \[
  -2A^{3/2}D_g Z=(\rho_g+p_g)+cy^3(\rho_f+p_f)\ge0.
  \]
  Interior root intervals have locally constant `y` and obey
  `D_g Z=-(rho_g+p_g)/(2G sqrt(A))<=0`. Remaining root points inherit
  the sign by continuity. No division by `P` at a root is used.
- If `P` is identically zero, use `D_g H_g=-(rho_g+p_g)/(2G)<=0`
  directly. In this case `Z` need not be monotone.

Thus, for any ordered times in the same regular interval,
`H_g(t1)<=0` implies `H_g(t2)<=0` for `t2>t1`; a strict initial negative
sign remains strict. This excludes contraction-to-expansion transitions
without an assumption on the order of vanishing of `H_g`.

The topology argument allows isolated roots, arbitrary branch switching,
accumulating root points and root sets with nonempty interior. It requires
the stated regularity, not a finite number of branch changes.

## Sharp CD mismatch

In physical `g` cosmic time, any parent in this class obeys

\[
\max\!\left(
 \left|H_{\rm parent}(-\tau/2)+\frac8{5\tau}\right|,
 \left|H_{\rm parent}(\tau/2)-\frac8{5\tau}\right|
\right)\ge\frac8{5\tau}.
\]

No `Hdot` error assumption is needed. A static regular Minkowski member
attains this error, so the endpoint threshold is sharp across the class.

## Evidence and exclusions

The [proof](notes/monotonicity.md) derives the arbitrary-lapse identity,
closed-root-set argument, endpoint bound and two exact distinction
controls. The [certificate](certificates/regular-hr-monotonic.json) replays
S6.4 and its lineage and checks independent Fraction coefficient identities.
The [source audit](notes/sources.md) gives the primary evidence boundary
without a global novelty claim.

No unweighted `H_g` monotonicity is claimed for interacting branches; an
exact local NEC background provides a counterexample. No `Z` monotonicity
is claimed for the zero polynomial; a decoupled de Sitter control rejects it.
Neither control supplies a vacuum or perturbation-health theorem.

Non-flat or non-common/non-bidiagonal geometries, singular or other
square-root branches, singular-endpoint continuation, quantum or
NEC-violating stress, nonminimal/derivative/shared matter, time-dependent
interaction coefficients and changed physical frames are outside the
theorem. General parent and UV matching, and P8, remain open.
