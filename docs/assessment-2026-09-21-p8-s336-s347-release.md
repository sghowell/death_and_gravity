# P8 release assessment: S336-S347 and the closure-driven matching work

Date: 2026-09-21. All acceptance runs for the twelve previously frozen
packets are now complete. **Original M/V/G/B/R/P8 remain OPEN.** Scoped
P8(a) and the frozen linear classification are unchanged.

This release publishes existing scoped results and the subsequent
MATCH-1/RATE4 research. It does not introduce another numbered scientific
checkpoint or infer closure from a larger test count.

## Scientific scope of the frozen packets

- S336 establishes an independent hard-radiative curvature direction
  invisible to the retained flat matching and soft data. S343 determines
  the first parity-even explicit-curvature contact basis and selected
  matter UV-pole cancellation. Neither fixes the finite parent coefficient.
- S337-S342 supply the selected local-tadpole, mixed-source, factorized
  bubble, quadratic-cancellation and triangle/box radiation calculations,
  including physical parameter contours and scoped finite remainder
  estimates. S340's zero is a selected real-TT quadratic cancellation,
  not a general statement that the full quantum radiation vanishes.
- S344-S347 determine known local degree-six curvature coefficients in
  explicit covariant lifts. Different lifts are converted before
  aggregation. The complete selected scalar-loop coefficient in S347 is
  negative, with its stated original normalized magnitude below 10^-207.
  This is a local coefficient of the classical limiting action, not a
  physical above-threshold approximation, a value for extra parent chi,
  a quantum decoupling theorem or a full curved/bounce matching result.

Each checkpoint's linked FORMULATION, proof notes and raw native report
remain authoritative for its assumptions. All 251 frozen source files and
all twelve reports are byte-for-byte unchanged. Their private, independent
fresh and repository source copies were compared again before release.

## Completed ordinary and CLI acceptance

All twelve checkpoints passed their separate original-SymPy ordinary and
CLI acceptance runs: **13,864 ordinary tests in total**. The runs used
three sequential batches, S336-S340, S341-S345 and S346-S347. Within each
mode, successors reused completed ancestor caches; their short times are
not independent cold-rebuild timings.

| Checkpoint | Ordinary tests passed | Pytest seconds | Separate CLI wrapper seconds |
|---|---:|---:|---:|
| [S336](../problems/P8/s6/continuation/s6_336/FORMULATION.md) | 710 | 13457.05 | 12329.530409708037 |
| [S337](../problems/P8/s6/continuation/s6_337/FORMULATION.md) | 796 | 13.96 | 14.441256207996048 |
| [S338](../problems/P8/s6/continuation/s6_338/FORMULATION.md) | 822 | 3.02 | 3.3729316249955446 |
| [S339](../problems/P8/s6/continuation/s6_339/FORMULATION.md) | 912 | 6.19 | 5.910190040944144 |
| [S340](../problems/P8/s6/continuation/s6_340/FORMULATION.md) | 933 | 0.84 | 0.5697968329768628 |
| [S341](../problems/P8/s6/continuation/s6_341/FORMULATION.md) | 1037 | 14430.33 | 12824.491466709063 |
| [S342](../problems/P8/s6/continuation/s6_342/FORMULATION.md) | 1006 | 10.64 | 10.6996804169612 |
| [S343](../problems/P8/s6/continuation/s6_343/FORMULATION.md) | 2124 | 5.46 | 5.326956582954153 |
| [S344](../problems/P8/s6/continuation/s6_344/FORMULATION.md) | 1265 | 2.99 | 2.449164832942188 |
| [S345](../problems/P8/s6/continuation/s6_345/FORMULATION.md) | 1566 | 12.66 | 13.242162834038027 |
| [S346](../problems/P8/s6/continuation/s6_346/FORMULATION.md) | 1355 | 19467.25 | 13821.436639582971 |
| [S347](../problems/P8/s6/continuation/s6_347/FORMULATION.md) | 1338 | 7.76 | 7.972734916955233 |

Every ordinary and CLI wrapper returned exit code 0 and checked that its
captured inputs and original SymPy implementation remained unchanged.
The full-precision ordinary-wrapper timings, input/batch hashes, native
report hashes and per-packet counts are in the
[machine-readable acceptance receipt](validation/p8-s336-s347-2026-09-21.json).

In particular, S346-S347's formerly pending acceptance is now complete:
2,693 ordinary tests and both CLI replays passed. The S346 ordinary run
finished at 2026-09-21 05:15:06 UTC, followed by S347 at 05:15:15 UTC.
This supersedes the earlier release hold without editing frozen packets.

## Covering full P8 snapshot

The full run including S347 passed **115,487 tests in 15,091.74 seconds**,
with all 945 captured test files present and unchanged. Its complete
source-preserving wrapper returned exit code 0 after 15,096.732481332961
seconds; all 7,012 captured inputs remained unchanged and the original
backend was restored.

Only this FULL run used the audited exact-GCD adapter. Its 128 independent
original-backend comparisons passed; final counters were 69,087 domain
fallbacks, 7,388 exact descents and 94 mixed fallbacks. Ordinary/CLI and
native packet construction retained original SymPy.

The decisive terminal receipts were:

```text
P8 collection complete: 115487 tests; all 945 captured files present and unchanged
115487 passed in 15091.74s (4:11:31)
S341_347_FULL_FINAL_EXIT 0 seconds 15096.732481332961 ALL_CAPTURED_INPUTS_UNCHANGED_ORIGINAL_BACKEND_RESTORED
```

Test snapshot SHA256:
`d6ef957e778af925761376433995f19b81af768608f5713be3a9b64405cf0539`.

All-input capture SHA256:
`03d163affef525c2eb18964a75b3c0769a7c2f32b5e268e7a63f13c000372b07`.

S341-S347 target-batch SHA256:
`cefaf250fd36b83b7f5dcdba6ab20a0c5c65a9936835a1be6356a3c98c045c3d`.

The terminal tool truncated part of the progress output. The collection,
final pytest result, adapter counters and complete-wrapper exit were
retained. No missing result is inferred from progress dots. The earlier
107,839-test snapshot ended at S342; it is not used to claim coverage of
S343-S347. Historical failed development attempts and the old 95,554-test
wrapper exit 1 remain failures, not relabelled successes.

## Distinct evidence level of the matching work

The [closure-driven plan](p8-closure-plan.md),
[matching decision audit](assessment-2026-09-20-p8-matching-decision-audit.md),
[finite-normalization obstruction](assessment-2026-09-20-p8-match1-renormalization-obstruction.md),
[parent-input comparison](assessment-2026-09-20-p8-match1-parent-input.md)
and [RATE4 v1 construction](assessment-2026-09-20-p8-rate4-candidate.md)
are included in this release.

Their three read-only diagnostics are run separately from the full frozen
P8 snapshot. They retain VERIFIED_N status for the stated algebra, interval
bounds and written assembly. The full regression does not promote them
to complete physical matching certificates.

RATE4's four new physical normalization conditions retain the known
matter-loop reference. Its first-order four-direction matching and
radiation-conversion budget are constructed, while the complete hard
remainder, complementary curved data and higher-order errors remain open.
The original parent is not rewritten or declared excluded.

## Preservation and next decision

The publication keeps each native hash and source manifest intact. A
path-scoped blank-at-eof whitespace exception is used only for immutable
packet files; new assessments, receipts, diagnostics and root documents
receive ordinary whitespace checks. No global setting is changed.
Unrelated P4/P9 work is excluded from the commit.

The next research milestone remains **RATE4.REMAINDER**: determine which
unfixed matching directions can enter the selected test, establish the
admissible analytic/infrared projection, and bound the complete subtracted
hard result and its errors against the lambda matching target. Curved
coefficients invisible or degenerate in the four elastic rates need their
own same-parent control before a bounce verdict. The existing plan's
conditional/inconclusive stop rules remain in force.
