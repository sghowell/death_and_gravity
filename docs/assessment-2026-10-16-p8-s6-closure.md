# P8 S6 closed under its own stop rule

Date: 2026-10-16. Decision: the user approved closing S6 after the
2026-10-16 portfolio review, which made P4 Theorem C the main track and P2
the one new problem. This note is a research-control record, not a new
certificate. No ledger row is changed or withdrawn. Original P8 remains
**OPEN**.

## 1. Final gate status

Copied unchanged from the closure ledger in
[p8-closure-plan.md](p8-closure-plan.md), as it stood at the last S6
commit (2026-09-29):

| Gate | Final status |
|---|---|
| A: scoped QEI/incompleteness objective | COMPLETE IN STATED SCOPE |
| L: frozen linear/matter classification | COMPLETE IN STATED SCOPE |
| M: decision-relevant quantum matching | OPEN; CONDITIONAL ONLY |
| V: vacuum necessary tests | OPEN |
| G: finite-gravity necessary tests | UNTESTED |
| B: same-parent bounce matching | OPEN |
| R: classification coverage and integration | OPEN |

Gates A and L belong to the earlier P8(a) and P8(b) work and are not
affected by this closure. Every gate that S6 owned (M, V, G, B, R) ends
open or untested.

## 2. Why S6 stops here

The plan's own **CONDITIONAL ONLY** stop rule says: stop expanding
known-sector calculations while the answer depends on unbounded physical
matching data, and resume only with a concrete new parent input or a
different route.

That rule fired three times:

1. **2026-09-20, MATCH-1.** A covariant order-hbar direction preserves the
   retained vacuum value and clock jets but changes the selected
   positivity combination C_pos
   ([audit](assessment-2026-09-20-p8-match1-renormalization-obstruction.md)).
2. **2026-09-21, RATE4.REMAINDER.** A local contact with at most ten
   derivatives leaves all four rate targets unchanged but shifts the
   analytic coefficient
   ([audit](assessment-2026-09-21-p8-rate4-remainder.md)).
3. **2026-09-24, PRESCRIPTION-2 and gravity.** Two explicit covariant
   functionals preserve the first-order matching data and reference mean
   while changing independent response directions. No quantitative
   same-observable gravity remainder exists
   ([assessment](assessment-2026-09-24-p8-response-closure-blockers.md)).

The S0-reference construction that followed on 2026-09-28/29 was a new
candidate. It was still IN PROGRESS and stated no decision threshold that
could move V, G or B.

S6 produced about 300 commits, about 350 assessment documents, and 365
ledger rows (S6 rows plus the MATCH-1, RATE4, PRESCRIPTION, COVZERO,
SOURCE and S0-reference studies). No top-level gate changed status after
2026-09-20.

## 3. What S6 established

- **S6.1 (vacuum-information independence).** An explicit smooth D-only
  family agrees with the original action on an open set containing the
  whole clock tube. It can carry either sign of the leading vacuum
  positivity diagnostic. Requiring global real analyticity instead forces
  the original F, which has no finite constant-clock Minkowski vacuum.
  Both are scoped algebraic statements
  ([S6 README](../problems/P8/s6/README.md)).
- **Many scoped intermediate results** (loop, soft, tree and matching
  computations) are recorded in CLAIMS.md with their own scope sentences.
  They remain valid as stated. None of them is a V, G or B verdict, and
  they should not be cited or counted as P8 results. Their CERTIFIED
  labels predate the WRITTEN_PROOF level, and none has had an independent
  review.

## 4. Status of the rest of P8

- **P8(b), linear classification (gate L):** complete in its frozen
  scope. A separate-agent adversarial review reported on 2026-10-16
  ([review](review-2026-10-16-p8b-adversarial.md)): the chain survives,
  with no claim-level error. Every load-bearing coefficient re-derives in
  an independent covariant expansion, and the replay passes. The gaps
  are presentation-level. The witnesses' contracting tails have p > 1/3,
  so homogeneous anisotropy grows there; the stated BKL/infrared
  exclusion covers this, but it belongs next to "stable".
- **S5, perturbative control:** parked. S5.1 is complete; there is no
  strong-coupling verdict.
- **S6, UV compatibility:** closed, as recorded here.
- **P8(a):** rows A.1–A.23. A.1–A.7 were reviewed in the Empiricist trial
  (unmerged branch `empiricist-claims-import`). Two separate-agent
  adversarial reviews covered the rest on 2026-10-16. Neither found a
  mathematical error, and every replay passed.
  - A.16–A.23 ([review](review-2026-10-16-p8a-adversarial.md)) are correct
    but narrow conditional FLRW theorems, built mostly on known QEIs.
  - A.8–A.15 ([review](review-2026-10-16-p8a-A8-A15-adversarial.md)) are
    internally correct but physically confined to Planck-scale durations,
    bounded by the classic SEE runaway, and do not advance P8(a).
  - Levels under the WRITTEN_PROOF policy adopted the same day: A.9–A.12
    and A.14–A.23 are WRITTEN_PROOF; A.8 (conditional on A.7's constants)
    and A.13 stay CERTIFIED.
  - The A.11 and A.15 test suites fail to load from long checkout paths:
    their conftest files put every `src` folder under problems/P8 on the
    import path. They pass from the main checkout.

The P8, S6 and P8(a) READMEs are hash-pinned by their replay
certificates, so this status is recorded here rather than in them.
Editing them would make `p8.verify_all --check` and the S6 replay fail.

### Validation receipts

The receipts in `docs/validation/` pin README.md, CLAIMS.md and
docs/p8-closure-plan.md by hash, as they stood at each receipt's own
checkpoint. Each receipt went stale at the next checkpoint. The 2026-10-16
ledger edits make the last one (2026-09-29) stale as well. They are
point-in-time records, not replays of the live ledger.

## 5. What would reopen S6

A reopening needs the user's approval, a stage budget, and a proposal
that meets the plan's admission rule. The proposal must name the open
gate, the decision threshold, why existing results cannot supply it, and
the next action for each outcome. It must also bring a concrete physical
input that fixes or bounds the finite matching combination entering the
chosen test, for example a specified UV parent with computed matching.
More known-sector refinement does not qualify.

## 6. Process note

Progress in S6 was reported by certificate and test counts rather than
gate changes for three weeks after the stop rule first fired. From now
on every stage carries a hard budget, and the coordinator checks gate
movement, not row counts.
