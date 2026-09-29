# P8 closure-driven work plan

Adopted 2026-09-20 following the user's request to make the next phase
explicitly closure driven. This is a research-control document, not a new
scientific certificate. Original P8 remains **OPEN**.

## Fixed finish line

The [original P8 statement](problems/open-problems-theoretical-cosmology-2026.tex)
and [adopted S6 contract](../problems/P8/s6/FORMULATION.md) remain authoritative.
This plan changes priorities, not their hypotheses or earlier evidence.

- Preserve the completed **free-photon/global-flat-FLRW P8(a) specialization**
  and its [explicit closure scope](../problems/P8/a/fields/maxwell/focusing/cosmology/notes/closure.md).
  Its nonminimal/null extensions remain separately qualified. Do not reopen
  it by demanding every field, arbitrary spacetime or an observed-universe fit.
- Preserve the completed [32-row linear-principal P8(b) classification](../problems/P8/README.md).
  Do not silently strengthen it to unrestricted nonlinear/BKL stability.
- Resolve the remaining vacuum, finite-gravity and common-parent matching
  obligations at the scope claimed. Necessary positivity is not sufficient
  for UV completion; constructing a complete all-orders UV theory is not an
  added completion requirement.
- A verdict for one candidate does not classify every allowed parent or
  operator row. Final integration must explicitly map the results to the
  M0 C-only, M0 D-only and M1 CD minimal-support branches and their supersets.
  A conditional or unresolved branch remains labeled as such. Do not rename
  a partial, conditional deliverable as closure of original P8.

## Closure ledger

| Gate | Present status | Evidence required to change status |
|---|---|---|
| A: scoped QEI/incompleteness objective | COMPLETE IN STATED SCOPE | Preserve the A.16-A.18 theorem, calibration and qualifications; no new prerequisite imposed here. |
| L: frozen linear/matter classification | COMPLETE IN STATED SCOPE | Preserve all 32 rows and their exact operator/matter hypotheses. |
| M: decision-relevant quantum matching | OPEN; original retained-data route and RATE4 remainder inference CONDITIONAL ONLY; RATE4's four-direction first-order normalization remains constructed | A physical common-parent determination or enclosure of the finite matching combinations that enter the chosen test, with a justified omitted-term budget and fixed physical frame. A subtraction convention alone is insufficient. |
| V: vacuum necessary tests | OPEN | A justified scalar observable/decoupling limit, pole/cut subtraction, applicable analyticity/unitarity input and bounded errors giving a definite necessary-test verdict. A formal first-loop representative alone is insufficient. |
| G: finite-gravity necessary tests | UNTESTED | The actual infrared prescription and a justified dispersion/contour error, combining the graviton pole and contour contribution before the forward limit. No guessed gravitational allowance. |
| B: same-parent bounce matching | OPEN | The same action, quantum prescription, state and matter frame must support the claimed bounce domain and completeness, with controlled heavy/omitted terms and relevant stability margins. A short finite regulated hybrid is not enough. |
| R: classification coverage and integration | OPEN | Map candidate/class-wide verdicts to every branch claimed in the final classification; distinguish existence from exclusion and necessary-test compatibility from UV completion. Independently review the full implication chain. |

M is a common prerequisite for a model-specific verdict, not a requirement to
compute every off-shell coefficient of an effective action. Derive only the
projections and bounds needed by V, G and B, unless a proof shows that the
larger calculation is necessary. These projections can differ; an on-shell
radiative result is not automatically a curved-background bound.

## First milestone: MATCH-1

**Question:** Can the unchanged QG2-H8A420 route supply physical bounds on
the finite coefficient combinations that already enter its proposed
positivity test, rather than only their known loop contributions?

The [first decision audit](assessment-2026-09-20-p8-matching-decision-audit.md)
finds that the existing data do not yet do this. Two explicit sensitivities
are nonzero; the separate S347 known-curvature calculation does not remove
them. This is a data-sufficiency failure, **not** an exclusion of the model.

The smallest first target is the local combination

    C_pos = alpha/(8*pi^2*kappa^2)
          + 4*g*(n+2*mu)*c_RH/[kappa*(n-2*mu)^3].

This is a projection of two already identified matching directions, not an
exhaustive formula for the physical positivity coefficient. The independent
Newton normalization, remaining sectors, radiative matching and omitted
orders must remain explicit wherever the chosen observable requires them.

### Current result and active subtarget

The [finite-normalization obstruction and four-value matching audit](assessment-2026-09-20-p8-match1-renormalization-obstruction.md)
exhibits a covariant order-hbar direction that preserves the retained
symmetric vacuum value and reference-clock jets but changes C_pos.
Even a four-jet clock-tube bound below 10^-590 permits opposite signs of
the selected tree-plus-matter-loop reference coefficient. This is not a
sign verdict for the full physical amplitude or a model exclusion.

The unchanged-retained-data inference route therefore receives
**CONDITIONAL ONLY**. Physical MATCH-1 remains OPEN. Its next subtarget,
**MATCH-1.PARENT-4**, is a specified quantum-parent determination of four
nonforward matched residuals at (s,t,u)/mu = (8,-2,-2), (10,-3,-3),
(12,-4,-4), (16,-6,-6), including the complete residual error. An exact
projector eliminates the unknown constant and Newton anchors and has
absolute error amplification below 15/mu^2 for n>=32mu. At the original
mu=1, epsilon<=lambda/15 would allocate at most lambda to this error;
no such physical input or error bound is presently supplied.

This algebra narrows the required calculation; it does not fill in missing
physical matching data. Compare and justify an additional parent input
before restarting numbered known-sector calculations. Keep the distinct
curved/state-dependent projection needed for B explicit.

The [parent-input comparison](assessment-2026-09-20-p8-match1-parent-input.md)
now supplies an exact conditional conversion from four finite-resolution
inclusive rates to C_pos. One explicit sufficient allocation bounds its
uncertainty by (243/320)*lambda, but the actual rates and complete error
bounds are still missing. At the previous amplitude precision, fitting
c_RH separately permits an error radius between 9*10^793 and 10^794;
the stable C_pos projection is not a curved-parent reconstruction.

The user has now approved developing the separately named,
observable-normalized **QG2-H8A420-RATE4** candidate. Its
[v1 specification and first-order construction](assessment-2026-09-20-p8-rate4-candidate.md)
adopt four rate targets equal to the retained matter-loop reference,
`R_i=1+2 Re L_m(q_i)/a_i`, at detector resolution `1/256`. These are new
physical normalization conditions, not old-parent predictions or measured
data. They preserve the matter reference rather than canceling it.

For fixed complementary inputs the four finite coordinates solve
`theta=-M^(-1)[Re F(q_i)+a_i D_i^rate/2]`. The rate Jacobian is invertible,
and the matched result is independent of finite hard-reference changes.
This completes the four-direction first-order normalization subtask, not
physical MATCHED status. The full retained Born-tree radiation conversion
is now bounded at first order: its contribution to the matching projection
is below `10^-783*lambda`. This does not bound the hard amplitude or higher
orders. The original one-loop hard-sector inventory and no-double-counting
rules are explicit; its complete subtracted bound remains unassembled.
Finite heavy-residue changes are exactly
degenerate with the curvature-exchange and constant directions in this
four-scalar amplitude; neither individual curved coefficients nor
complementary higher-EFT terms are fixed by the four rates alone.

The [RATE4.REMAINDER decision audit](assessment-2026-09-21-p8-rate4-remainder.md)
now gives **CONDITIONAL ONLY** for this candidate family's complete
projection bound. A covariant local contact with at most ten derivatives
leaves all four first-order rates and fitted coordinates unchanged but
shifts the analytic coefficient by `-902256*zeta`. The missing bound is
decision-relevant even with exact known loops. The audit also bounds the
entire gapped H/Proca subtracted increment below `10^-1003*lambda` and
separates the heavy-residue cancellation from mass-shift sensitivity.
This does not complete the full hard bound or justify its massless V/G
functional; the original parent and RATE4 v1 remain conditional.

The [RATE4.PARENT comparison](assessment-2026-09-21-p8-rate4-parent.md)
has now examined the additional-input routes. It recommends a separately
named **RATE4-D8, parent-bounded** proposal retaining the approved targets.
Its complete eight-derivative local fit has coefficient error norms below
27 and a sufficient conditional total error below
`(1620429/3259000)*lambda < lambda/2`. If the stronger tail restriction on
SU is inappropriate, **RATE5-D10** adds one off-angle rate and fits that
coefficient; its local error norms are below 20 and its conditional total
error is below `(987/2000)*lambda < lambda/2`. Both require independent
pole/source and weighted-tail bounds. These are proposed additional inputs,
not established properties or adopted conditions of RATE4 v1.

The eight-rate alternative is algebraically invertible and has a stable
target projection, but at independent sample error `lambda/100` its J
coordinate radius lies between `10^1168` and `10^1169`. Stable b20
inference must not be confused with controlled individual heavy or curved
matching. This is conditioning, not an exclusion of exact eight-rate data.

The user has now approved developing the separately named
[RATE4-D8 v1 conditional candidate](assessment-2026-09-21-p8-rate4-d8-candidate.md).
Its retained four targets, complete-input matching and interval maps,
conditional formal-order normalization lemma and sufficient analytic-tail
witness are specified. Approval authorizes this conditional definition,
not the existence of a physical parent satisfying its bounds. RATE5-D10
and its fifth target remain unadopted; no original frozen calculation changes.

**RATE4-D8.WITNESS result: retained-source route insufficient.** The
[source and identifiability audit](assessment-2026-09-22-p8-rate4-d8-witness.md)
finds no actual parent-bound witness. Entire local S^k R5 comparisons
preserve any finite weighted S/U coefficient prefix, the four first-order
rates and known first-order cuts while leaving the unknown tail and its
forward projection unbounded. Separate Q/J/T compensations vary pole/source
coordinates without changing the analytic tail. These are formal EFT
data-sufficiency obstructions, not physical UV counterexamples or exclusions
of a parent that actually satisfies D8. Stop the unchanged-source inference
route; another finite coefficient prefix cannot close this gap by itself.

**PRESCRIPTION-1 feasibility study complete:** the user authorized all three
avenues: a concrete quantum prescription, a smaller decision-projection
bound and gravity/bounce compatibility. The
[result and written proofs](assessment-2026-09-22-p8-prescription1-feasibility.md)
give scoped negative and positive decisions, not a complete quantum parent:

- A finite proper-time covariance fails a necessary positive-spectral and
  reflection-positivity test when declared to be the exact physical local
  scalar covariance. A sharply band-limited field is not the original exact
  local field. These identifications stop; computational regulators and
  Wilsonian EFT remain possible.
- The selected light-bubble omitted short-time projection is bounded by
  `33*abs(beta)/Lambda^10`. A restricted positive scalar-exchange class is
  controlled by one physical spectral moment. Neither determines the
  independent local/gravity boundary action. RG compensation explicitly
  preserves that free physical integration constant.
- The actual S240 heavy-state renormalized readout tail is below `10^-697`
  above a test fixed comoving split `P=10^100`, uniformly through five time
  derivatives on `[-1,1]`. The state and all local terms are retained. This
  is not an interacting response, physical mode cutoff or global B theorem.
- The matching projector is not a positive gravity functional: it has a
  negative value even on the tested positive scalar exchange. No admissible
  finite-gravity error allowance or common interacting bounce is supplied.

**Next admission boundary:** a proposed local, source-bounded matching model
must determine or bound the remaining independent local/gravity and curved
terms, not merely add another negligible heavy sector. Action-determined
positive-spectrum correlations are one surviving mechanism; the study's
scalar example does not complete the original parent. A fully specified
fixed-scale quantum-EFT boundary definition can also be proposed, but its
finite values are new physical assumptions, not subtraction conventions or
predictions of the present classical action. State their effect on the four
targets, operator support and prepared state before adoption.

Keep the current H8A420/RATE4-D8 route conditional and parked for full
matching until such input is supplied. Do not repeat the already authorized
feasibility request or the excluded regulator identifications. An actual
source-bound model can resume M, followed by the independent V/G and B
tests; a partial spectator bound alone cannot. This is a scoped research
boundary, not an all-parent no-go or a new all-orders UV-completion
requirement. A theorem for the actual finite-gravity observable remains an
alternative only if its hypotheses and absolute error are verified.

**PRESCRIPTION-2 and gravity applicability are now started.** The
[first constructive packet](assessment-2026-09-22-p8-prescription2-gravity.md)
contains a proposed RATE4-COVZERO-1 additional finite-functional rule,
explicitly unadopted physical boundary data rather than a scheme inference.
Its clock-localized lift has no reference background or quadratic-response
change at first quantum order. This is a difference-functional jet lemma,
not a complete interacting B theorem. A common covariant renormalized
reference still has to be constructed; flat component references cannot
be glued together by assertion. A local heavy-source alternative supplies
flat/curved coefficient correlations but does not determine independent
gravity matching by itself.

In parallel, an actual ball-autocorrelation weight now gives a massive-arc
positive comparison for every partial-wave spin and s>=32. At unchanged
detector1/256, the finite-gap and crossed-denominator correction is below
J/80, with J a separately required positive forward spectral moment. A
finite bound J<=20lambda would allocate lambda/4 to this contribution;
it is a sufficient target, not an adopted parent condition. The original
leading heavy-scalar part has J_H<5lambda, not a full J bound. Complete
finite-coupling/unitarity, contour and low-cut errors remain separate.

The arc differs from the forward functional. Its freshly derived RATE4
weights have norm below one and its D8 tail error is below lambda/400 IF
the old tail hypothesis holds. Exact finite-gap positivity fails even at
spin256 in the displayed control; do not replace the all-spin proof and
error estimate with a finite low-spin search. The sixth matching moment
alone does not control the required third forward moment.

**September 24 continuation:** the
[reference-reduction and finite-window packet](assessment-2026-09-24-p8-reference-regge.md)
isolates the 56 algebraic connection modes and their ultralocal factor
in a stated dimensional component-measure convention; it does not supply
the constrained physical measure. The original heavy source has no direct
formal one-loop reference-mean or quadratic-response contribution, but
its free metric determinant remains. The original Proca source has a
nonzero quadratic jet and produces a quartic exchange with a nonzero
spatial-momentum control. That response cannot be dropped using a
homogeneous or source-free Gaussian calculation.

On gravity, an explicitly positive even-spin spectral fixture has
divergent full J but a finite finite-angle arc. The comparison can instead
be restricted to a finite energy window, with a separate finite-angle
Regge envelope bounding the tail by <60 C/alpha plus its supplied
remainder. This does not establish that the actual parent's J diverges
or supply its actual Regge/moment/contour/observable inputs.

**Further September 24 progress:** the
[constrained-source/pole-matching packet](assessment-2026-09-24-p8-source-response-pole-match.md)
uses the already existing S250 filtration, S253/S257 finite-regulator
constraint measure and S275 boundary correction. These are inherited
results, not newly solved continuum-reference obligations. The actual
time-dependent source cubic and the lapse-generated source-squared
quartic are now derived in the original prepared coordinates. Their
scalar and longitudinal-vector one-loop kernels include both contacts
and distinct momenta on each leg. This is a formal tree-clock-vertex,
same-state finite-regulator result, not an evaluated continuum response.
The contact is indefinite and cannot be silently omitted or normal-ordered
to a new physical zero.

On gravity, single-leading-pole matching in an explicitly restricted
dispersion scope fixes C(0)/j'(0)=pi/(2kappa). It does not supply the
finite-angle residue/trajectory shape bounds needed to use that value
in the tail estimate. The finite part also depends on residue and
trajectory derivatives. Published massless finite-energy sum rules do
not determine these data from the four RATE4 samples or close the
actual observable's massless-loop/contour problem.

**Closure-driven September 24 outcome: physical-input stop.** The
[continuum-extension and closure audit](assessment-2026-09-24-p8-response-closure-blockers.md)
establishes the same-state source-squared noise distributions and
time-retarded extension families. The conservative complete-phase
scaling degree is <=14; local normal-derivative ambiguities have order
at most10. Existence is not a chosen covariant subtraction or B bound.
Two explicit regular covariant finite functionals preserve the retained
first-order flat normalization data and reference mean while changing
independent reduced response directions. Their response determinant is
102400/59049. They are unadopted finite quantum-completion witnesses,
not modifications of frozen sources or all-orders equivalent theories.

The current gravity route also cannot certify its numerical allowance
from the retained inputs. An O(G) finite-detector estimate has no supplied
constant or validity domain for this exact observable. In the original
normalization an allocated error C_err/kappa<lambda/100 would require
C_err<10^198 and applicability at the actual parameters. Neither that
physical bound nor the finite-angle/contour and observable conversions
has been proved. No Regge factor or missing error is assigned a value.

**Decision resolved by user approval:** RATE4-COVZERO-1's first-order
finite-boundary rule is now adopted as a separately named candidate, not
a consequence of the four rates. The
[reference-seed packet](assessment-2026-09-24-p8-covzero-reference-seed.md)
constructs a dimension-dependent light/metric action with zero kinetic
Schur complement at arbitrary tilt, exact clock ADM cancellations and
a regular source/chart through X=0. A fixed four-dimensional coefficient
at all d fails this constraint-preserving check. The new continuation's
evanescent terms must be kept through pole subtraction, not assigned zero.
The explicit Proca mu=1000 to mu=1 local conversion preserves the old
finite prescription, which a scale change alone would alter. Auxiliary-
gradient cotangent bookkeeping addresses the chart, not the full quantum
gauge measure. The candidate is recorded in a machine-readable
[specification](candidates/p8-rate4-covzero-1.json).

**Reference admission result: exact fixed-count extension obstructed.**
The [coupled-source audit](assessment-2026-09-24-p8-covzero-source-constraint-obstruction.md)
keeps the full retained source and finds a generic scalar/metric rank
increase from six to seven. At a purely spatial scalar gradient the
new acceleration Hessian is (R-1)^2/X, and the temporal-Proca secondary
equation still solves its own auxiliary variable instead of replacing
the lost scalar constraint. A literal-parent bound for u=0,
X=-e,0<e<=10^-430 gives 5*1024*e^3<R-1<10*1024*e^3 and a nonzero
metric pivot. This is an off-shell domain obstruction to the proposed
uniform exact canonical measure, not a demonstrated bounce instability
or a no-go for every controlled EFT. The b=0 clock count, vacuum
quadratic spectrum and pure-light reference-seed results are unchanged.

**SOURCE-1: authorized bounded study COMPLETE; subsequent S0 adoption approved.**
The [source feasibility study](assessment-2026-09-28-p8-source1-feasibility.md)
independently finds a nonzero scalar/vector primary bracket at nonzero null
gradients in the actual parent. Within the fixed-light, Hessian-linear
gradient-parallel source class, the two-tilt constraint test forces the
second-derivative coefficients to vanish on the regular set. Source-free
S0 and a first-derivative-only comparison have safe primary algebra and regular
local secondary pairs at the vacuum and finite clock points. They preserve
the classical clock and quadratic data. A field-degree proof preserves the
formal source-tagged first-loop flat four-scalar projection, but the curved
nonlinear response changes. The unchanged-action order-reduced EFT route
still needs an actual validity/error bound; it is not excluded.

**Decision resolved: S0 is now adopted and S0-REFERENCE is in progress.**
The user explicitly approved the separately named
[QG2-H8A420-RATE4-COVZERO-S0 candidate](candidates/p8-rate4-covzero-s0.json),
replacing the entire source square by the source-free Proca mass term.
The old sourced parent is preserved. The
[canonical/vector bridge](assessment-2026-09-28-p8-s0-reference.md)
establishes compact-clock local auxiliary admission, the exact source-free
four-plus-two Gaussian split and reduced canonical brackets, physical
vector probe contacts, original-state transport and the inherited ordinary
Proca finite conversion. The old source-aligned response difference contains
both the two-vertex kernel and its local contact; deleting the source does not
delete the vector's gravitational response. A restricted lapse-only
Horndeski-frame shortcut is singular at finite clock times, while the
original canonical chart stays regular.

The [gauge/measure continuation](assessment-2026-09-29-p8-s0-gauge-measure.md)
now gives the dimensional clock-gauge algebra and temporal Jacobian,
the actual de Donder ghost operator, its momentum-uniform finite-interval
retarded bound and ordinary ghost UV polynomial. It also corrects the
previous unqualified held-vector response comparison: retain the old
normal-vector one-point embedding contact `-Gamma_old,W*S_AB`. A centered
vector product does not remove its light quadratic contribution. This
correction does not change the adopted S0 action or its parity. Original
packets stay frozen; successful replay is not endorsement of that
superseded interpretation.

**Active subtarget:** construct the coupled S0 light/metric/second-class
measure bridge in the selected dimensional continuation, with actual
off-shell source/state and endpoint transport. Combine its bosonic and
constraint contributions with the ghost block; establish the needed
Ward/locality conditions and common pole/finite conversions. On-clock
determinant agreement alone is not the response bridge. Then
assemble the same-state causal finite parts and actual response bounds.
Compact local/patchwise admission is sufficient if it meets the original
closure contract; no unnecessarily global field-space prerequisite is
added. The scalar-unitary chart is not extended through the flat vacuum.
No new physical domain restriction, state, cutoff or mode prescription is
adopted. Both S0 and COVZERO approvals are resolved; do not request them
again or ask the user to invent Wilson/Regge coefficients. The common
reference remains incomplete, not an executable certified quantum parent.

After resolving this admission boundary, finish the common reference and
evaluate/bound the response of the named candidate. For the unchanged
source this includes its source-squared kernels, local contacts and mixed
response; for S0 those source-tagged vertices are absent, but all retained
source-independent light, metric and matter response remains. Independently
establish a same-observable gravity certificate with quantitative errors
and a valid domain; authorization is not proof of those hypotheses.
The full-J route remains optional, and a suitable finite-energy theorem
is admissible. The earlier inference routes from retained data alone
remain CONDITIONAL ONLY; the authorized COVZERO candidate is additional
physical data, not a proof of the missing gravity bounds. Stop expanding
known-sector refinements as a substitute. This is not exclusion of all gravity routes
or closure of M/V/G/B/R/P8. No background research is claimed after handoff.

The optional R=128 bidisk witness (proved supremum <=23000lambda) remains
sufficient, not necessary. A direct tail or decision-projection bound may
be preferable. Pole/source enclosures must come from the same parent's
physical prescription. Retain the distinct curved/source and original
prepared-state obligations for B, the complete hard calculation, actual
higher-order errors and admissible V/G functional. The order-by-order
normalization lemma is not convergence or an omitted-order bound. The only
newly adopted physical data are the authorized COVZERO first-order boundary
rule and the separately named S0 source replacement; no physical cutoff or
weakened finish line is adopted. Keep
CONDITIONAL ONLY until actual bounds are established; do not
ask the user to invent Wilson coefficients or restart negligible known-
sector refinements as a substitute.

### Required output

1. Identify the parent/prescription input that fixes or bounds this
   combination. Separate coefficients fixed by a physical condition from
   those whose finite value is merely a convention or remains unknown.
2. Derive a finite enclosure for C_pos in the same normalization as the
   chosen observable, or prove a precise obstruction for the specified
   matching ansatz. Do not infer an enclosure from small known loops.
3. State the remaining error terms and the proof needed to bound each one.
   Give the decision threshold they must satisfy; do not report a small
   number without its comparison margin.
4. Show where the same finite data enter B, or explicitly identify a
   different required curved projection. Retain source/matter contacts in
   any field redefinition.

### Outcomes and stop rules

- **MATCHED:** a physically justified enclosure is obtained. Proceed to a
  bounded V/G test and B compatibility using that same input. This is not
  itself a positivity pass or P8 closure.
- **EXCLUDED ANSATZ:** a proved contradiction rules out this specified
  matching route. Record the exact hypotheses and pivot; do not infer an
  all-parent or all-row no-go.
- **CONDITIONAL ONLY:** the answer still depends on unbounded physical
  matching data. Stop expanding known-sector calculations as if they could
  resolve that ambiguity. Identify a concrete new parent input or compare
  a different route before resuming a numbered sequence.
- An enclosure straddling a decision threshold is **INCONCLUSIVE**, not a
  failed physics test. Refine only if there is a specific bound likely to
  change that verdict.

No arbitrary alpha, c_RH, chi or Regge constant is to be requested from the
user to manufacture a result. If a separately named candidate introduces
physical conditions, label them as new assumptions and test its full
matching implications rather than rewriting the old candidate.

## Sequence after MATCH-1

1. **V/G decision:** select one necessary inequality and its actual physical
   observable. Assemble a signed coefficient enclosure and a justified
   gravitational allowance; retain all relevant truncation/infrared errors.
   A finite-energy contour theorem is sufficient if its hypotheses and
   absolute error can be established. An asymptotic all-energy theorem is
   not mandatory solely because earlier notes listed it as future work.
2. **B compatibility:** establish that the *same* matched data admit the
   matter-coupled bounce with the required domains, tails and error margins,
   or give a scoped obstruction. Do not substitute a different state,
   frame, finite prescription or independently healthy vacuum.
3. **R coverage:** propagate only justified results through the operator-row
   classification. Identify any uncovered C-only/D-only/CD branch before
   writing a closure conclusion.
4. **Final review and release:** replay relevant evidence, check that no
   unresolved assumption carries the conclusion, and update the claims
   ledger and public status together. No status is inferred from test count.

## Admission rule for further calculations

Before starting a new scientific checkpoint, record:

- the specific open gate and missing term it addresses;
- the decision threshold or obstruction it could establish;
- why existing results cannot supply it;
- the next action for success, exclusion and inconclusive outcomes.

Do not add another certificate merely to tighten an already negligible
known contribution while the decision is limited by unbounded matching or
contour data. Keep genuine partial results, but do not count them as closed
top-level gates. Report progress by gate changes and remaining assumptions.

## Separate validation/publication lane

The [2026-09-21 release record](assessment-2026-09-21-p8-s336-s347-release.md)
now records completed acceptance for S336-S347, alongside this plan and
the MATCH-1/RATE4 diagnostics. Their scientific files and stored
certificates remain unchanged. Publication acceptance is not M/V/G/B closure.

- Retrieved the successful 107,839-test full regression on 2026-09-20.
  It ran through the snapshot containing S342 and does **not** cover
  S343-S347. Runtime was 15,333.71 seconds; the source-preserving wrapper
  returned exit code 0.
- S336-S340 ordinary/CLI acceptance was already complete. S341-S345
  ordinary acceptance (6,998 tests total) and all five separate CLI replays
  are now confirmed complete with exit code 0.
- S346-S347 ordinary/CLI acceptance has now passed: 2,693 ordinary tests
  and both separate CLI replays, with complete-wrapper exit code 0.
- The full snapshot including S347 passed **115,487 tests in 15,091.74s**,
  covering 945 captured test files. Its wrapper returned exit code 0,
  preserved all 7,012 captured inputs and restored the original backend.
  The audited exact-GCD adapter was FULL-only; ordinary/CLI used original
  SymPy. The three root matching diagnostics are verified separately.

The validation wait is resolved. The subsequent RATE4.REMAINDER audit is
a separate read-only diagnostic, not another full-suite replay. Proceed
with PRESCRIPTION-2 reference completion and the gravity spectral/error
calculation after their initial constructive packet, preserving unrelated
P4/P9 changes. Never
modify a frozen certificate to make a replay pass or use release counts
as a substitute for closing a scientific gate.
