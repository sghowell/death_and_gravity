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
| M: decision-relevant quantum matching | OPEN; original retained-data route CONDITIONAL ONLY; approved RATE4 candidate has a constructed four-direction first-order normalization | A physical common-parent determination or enclosure of the finite matching combinations that enter the chosen test, with a justified omitted-term budget and fixed physical frame. A subtraction convention alone is insufficient. |
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

The active calculation is **RATE4.REMAINDER**: assemble and bound the
complete physically subtracted projection, with an admissible V/G
functional and the earlier absolute `lambda` matching-error target.
Separate calculated terms, the four removed directions, and independent
complementary data; retain the distinct curved/source projection for B.
No new arbitrary coefficient is requested from the user, and the original
parent stays conditional. Do not resume numbered known-sector refinement
as a substitute for this full remainder bound.

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

The validation wait is resolved. Proceed with RATE4.REMAINDER according
to the admission rule, while preserving unrelated P4/P9 changes. Never
modify a frozen certificate to make a replay pass or use release counts
as a substitute for closing a scientific gate.
