# P8 closure audit after A.10 and S6.5.HR

Recorded 2026-09-06. This is a status and research assessment, **not a new
theorem, certificate or amendment of the original problem**. P8 remains open.
The original source is
[`open-problems-theoretical-cosmology-2026.tex`, P8](problems/open-problems-theoretical-cosmology-2026.tex).
The frozen 32-row linear-principal classification is complete; that scoped
entry-point success does not close the realistic-field and UV questions.

## New completed gates in this continuation

| Gate | What is established | What is not implied |
|---|---|---|
| [A.8](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/FORMULATION.md) | Actual prepared-state density, pressure and EED, including anomaly and both retarded histories; quantitative actual SEE residual | A nearby exact semiclassical solution |
| [A.9](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/qsei/FORMULATION.md) | Every Hadamard target state and every compact real proper-time sampler in the specified region obey the bound with coefficient `5*hbar/(16*pi^2)` | Transfer to an unknown solution, all relevant geodesics or realistic field content |
| [A.10](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/FORMULATION.md) | Actual nonlinear mode continuity in the first-potential-derivative norm, common-history shortening, and an isolated causal logarithmic inverse | A complete self-mapping SEE contraction, compatible initial constraint, smoothness or interval shadowing |
| [S6.3.beta1](../problems/P8/s6/matching/bimetric/FORMULATION.md) | Positive quadratic vacuum spectrum and source-aware flat tree matching for the specified parent; exact CD obstruction | Nonlinear/rolling matching or arbitrary-parameter vacuum health |
| [S6.4.HR](../problems/P8/s6/matching/bimetric/general/FORMULATION.md) | Arbitrary-parameter, root-covered nondegenerate flat-bounce obstruction with separate NEC matter | Other matter couplings, geometries or general UV exclusion |
| [S6.5.HR](../problems/P8/s6/matching/bimetric/general/monotonic/FORMULATION.md) | All-degeneracy contraction-to-expansion exclusion in that class; sharp endpoint-only CD mismatch | A rejection of all parent mechanisms or an accepted UV matching construction |

The latest input hashes are:

    A.10: c37ee6cb74eb8383abbddf2f240c28676671e57b9f114dd8ec2195267de2d3ce
    S6.5.HR: 0efd16fa4f02f45f35448e2056fb37c863c29c2ed09203c3351ba9515baff031

Both are children of unchanged replayable certificate chains. Their own
formulations, source hashes, exact arithmetic and independent audits delimit
the claims; this assessment must not enlarge them.

Validation at this checkpoint: all **1,096 P8 tests passed** in the
joint run, including 74 A.10 tests and 23 S6.5.HR tests. Both new
read-only certificate replays and the full P8 Ruff check pass. The earlier
999-test checkpoint was independently committed as `6618b95`; A.8/A.9
and the beta1 source screen were committed as `f65e911`. Unrelated P4/P9
working-tree files were not included in these commits.

## P8(a): the next decisive task is a full applicability argument

The main missing result is no longer a coincident-stress calculation or
the nonlinear mode functional's derivative continuity. The chosen route
now needs a complete actual semiclassical solution construction or
validated shadowing argument. Its concrete completion conditions are:

1. A specified positive-metric function space with compatible actual state,
   gravitational initial data and classical radiation constraint.
2. The full integrated-trace map, retaining its geometric, anomaly,
   renormalization and gravitational terms. The A.10 nonlinear coefficient
   is only one term in its Lipschitz bound.
3. A genuine self-map and strict contraction or another existence estimate,
   followed by an enclosure over the entire intended target interval.
   An arbitrarily short local estimate is not a cosmological-duration bound.
4. Smoothness and Hadamard compatibility, followed by a two-point-function
   and reference-stress estimate on the resulting metric and geodesic domain.
5. Verification of the focusing theorem's initial geometric/domain hypotheses
   and a field-content calibration that addresses the original realistic-field
   requirement rather than only one massless minimal scalar.

The causal inverse has a positive exponential pole contribution; it cannot
be omitted to manufacture long-time stability. Similarly, freezing the
old metric and state up to an initial slice fixes its initial constraint
value. An off-shell mismatch cannot be corrected by declaring a new energy
value or by changing metric jets without tracking the state.

The [Meda--Pinamonti--Siemssen existence paper](https://arxiv.org/pdf/2007.14665),
Theorem 5.9 and Remark 5.2, supplies a useful local mild-solution framework
with explicit compatibility assumptions. It is not a ready-made smooth
massless shadowing theorem for our prepared metric. The
[Gottschalk--Siemssen framework](https://arxiv.org/pdf/1809.03812v3),
Theorems 5.1, 5.6 and 5.11, provides moment-space evolution, continuous
dependence and a compatible smooth preparation construction under its own
hypotheses. Applying it here would still require a justified realization
of the actual transported state in the required space and a quantitative
comparison to the chosen geometry. Setting formal initial moments to zero
is not such a realization proof.

These are alternative tools, not contradictory existence claims. A.10
avoids a mass gap for its own actual-mode estimate; that does not delete a
frequency or state assumption from a theorem proved using another method.

An exact solution near A.8 is the currently selected applicability route,
not a logical prerequisite for every possible new incompleteness theorem.
A different rigorous realistic-field/geometric theorem could address the
original P8(a) directly. It would need its own nonvacuous hypotheses and
calibrated constants; it must not be represented as an automatic consequence
of the fixed testbed.

### Unpromoted next-map algebra

A direct exploratory substitution into A.8's unchanged
`reconstruction.trace` gives a useful candidate organization. In actual
conformal time let `h=a'/a`, `U=-a''/a`,
`S=Rmode-K/2+U/10` and `Q=S/a^2`. With the named zero finite-curvature
coefficient and traceless classical radiation, the traced equation
`R= kappa*trace(T)` can be written

\[
 Q''+2hQ'=\frac{48\pi^2}{\kappa\hbar}U
  -\frac1{a^2}\left[\frac{U^2}4+\frac{h^2(U+h^2)}{30}\right].
\]

The direct symbolic difference from the pinned actual trace was zero
in this review. This is recorded as next-map algebra, not a new
certificate or solution theorem. With `a=A*eta_star*atilde`,
`x=eta/eta_star` and physical `epsilon=1`, the dimensionless coefficient
of `U_dimensionless` becomes `1/(60*delta)`, using A.8's coupling identity
`kappa*hbar=2880*pi^2*delta*A^2*eta_star^4`.

Thus A.10's small nonlinear-mode coefficient does not control even
this elementary term of the full map. The large coefficient alone
does not prove instability or a necessary maximum step length; a valid
inverse, cancellation or branch-selection argument must treat it explicitly.
The auxiliary initial values also have to equal those of the actual state,
and the initial density constraint remains a separate equation.

## P8(b): stop treating excluded parent families as general UV closure

The [frozen classification](../problems/P8/FORMULATION.md) has already decided
all 32 M0/M1 rows, including all-time physical completeness and the stated
linear-principal health and speed conditions. The later finite-band,
finite-order interaction and local-operator bounds are useful extensions,
but they do not by themselves prove every stronger meaning of stability.
Any claim of nonlinear or all-wavelength stability must identify that
stronger contract and supply the corresponding estimates.

The [adopted S6 framework](../problems/P8/s6/FORMULATION.md) is still missing
controlled matching between an admissible vacuum/parent description and
the surviving physical bounce, with the relevant extra EFT operators and
finite gravitational/IR allowances accounted for. The bimetric results
are now strong enough that another search over constant HR parameters
with the same regular geometry and separate NEC matter cannot rescue
that route, even by using a degenerate bounce or shifting the bounce time.

They do not decide derivative/nonminimal/shared matter, other geometries,
quantum stress or general parent constructions. A next parent proposal
must identify exactly which excluded hypothesis it changes, derive its
full physical-frame equations before truncation, and test background
compatibility before investing in a flat-vacuum expansion. If it survives,
its gap, retained modes, matching remainder, matter coupling and necessary
positivity conditions need quantitative control on the actual bounce.

A smooth spliced healthy vacuum is not controlled matching. Conversely,
the original question about consistency with UV positivity must not be
silently replaced by a demand for construction of a complete all-orders
UV theory: necessary positivity with justified applicability is distinct
from a full UV completion. Either a valid surviving realization or a
no-go covering the specified allowed class would be decisive; the finite
list of excluded parent families supplies neither general verdict.

## Closure decision

The six new gates are completed as scoped. Original P8(a), the outstanding
S6/positivity applicability work in P8(b), and P8 overall are **not closed**.
Closing the project now would either drop original obligations or promote
conditional/partial results to claims they do not prove.

No missing credential, approval, paid compute allocation or user preference
has been identified. The outstanding obstacles are mathematical. They
must be resolved by further proof and validation, not by changing a status
label or repeatedly asking the user to approve the same research direction.
