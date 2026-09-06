# Level semantics of this ledger (death_and_gravity)

From this repository's README (the level table in force since the first certificate):

| level | meaning |
|---|---|
| HEURISTIC | numerical/analytic lore, no certificate |
| CONJECTURED | precise statement, tested numerically, not certified |
| VERIFIED_N | independently re-run by a second implementation/engine |
| CERTIFIED | machine-checkable certificate (interval enclosure, exact dual, exhaustive enumeration) committed under `certificates/` |
| FORMALIZED | kernel-checked proof (Lean 4) |
| REFUTED | terminal |

As kept by the Empiricist claim ledger on this branch:

- A level names the strongest evidence on record for the exact statement text; it is
  not a truth value. `level` rises only through `promote`, which replays the evidence
  with a certified verifier (the `<package>.verify` checkers through
  `tools/empiricist_check.py`) and requires PASS.
- CERTIFIED here means: the committed certificate is reproduced exactly by the checker's
  read-only replay (`build_report()` compared with `validate_report`), AND an independent
  review receipt with no blocking finding attests that the certificate establishes THIS
  statement within the verification boundary the certificate discloses. That boundary
  is disclosed differently across certificate schemas: a `verification_boundary` field
  where present (later P8(a) certificates), otherwise the `status`, `not_established`
  and `written_proof` fields together with the named notes file. Written analytic steps
  the certificate places outside machine replay are permitted at CERTIFIED provided they
  are disclosed in one of those places; they are not machine-checked. Warnings in a
  receipt (a REVISE verdict) do not veto CERTIFIED; they stay on record for the author.
  FORMALIZED is reserved for kernel-checked proofs.
- A claim's level may not exceed the lowest level among the claims it depends on
  (`depends_on`); pinned prior certificates are locked path dependencies, so a change to
  a prior certificate makes the claim STALE.
- A level imported from the pre-ledger table is kept as `legacy_level` ("not
  re-earned") until `promote` re-earns it; `IMPORTED` evidence entries record which
  files that table cited and never count as verification.
