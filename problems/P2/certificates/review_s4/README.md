# review_s4: the independent review's checks of P2 S4 (part (i) at P_B)

These are the scripts and outputs behind `docs/review-2026-10-18-p2-s4-adversarial.md`.

- The review was run at `0e329217`, from a scratch directory, on 2026-10-18/19.
- The follow-up check of the author's fixes was run at `21358225` on 2026-10-19. Its outputs are in
  `followup_2026-10-19/`.
- The scripts were committed here after the review, with the user's approval. Their absolute paths
  were replaced by paths relative to the script, so they replay from any checkout.

**Running them.**
- Run each script from this directory with the repository's Python environment (`.venv/bin/python`
  at the root of the main checkout): `python <script>.py <out>.json`.
- Every script writes the JSON file named on its command line.
- The `.txt` files are the console logs of the original runs (they were `.log` files; the repository
  ignores `*.log`).
- Runtimes are single-process CPU times on the review machine.
- "No `p2`" means the script imports nothing from `problems/P2/src/p2`.

**Re-running at a later commit.** The outputs of the first review (top level) were computed against
the certificates as of `0e329217`. Two kinds of script compare against the stored certificates:
- the replays (`replay_author.py`, `replay_all_SL.py`);
- the independent checks that read stored bounds (`indep_S.py`, `indep_L.py`, `indep_TA.py`,
  `cover_recheck.py`).

Re-run at a later commit, they check the certificates of that commit. That is what
`followup_2026-10-19/` holds.

## The independent checks (no `p2`)

| script | what it checks | result (at `0e329217`) | runtime |
|---|---|---|---|
| `indep_params.py` | P_B constants from the exact definitions (mpmath, 50 digits); imported by the others | κ₋, κ₊, κ_c and γ₀/κ₋ agree with P2-par; X₀ = −2.2966368932; a* = 0.0299686149 | < 1 s |
| `indep_TA.py` | r(s) on S2's contour (scipy DOP853); for every T box, the bound E₀/h₀ + E₁ + h₀ min(0, E₂) on a 1/256 grid with tails; S2's E from V_l directly at actual l; A, and its monotonicity in Re ω | every certified μ ≤ the float minimum; hE ≥ μ_cert at 54 samples (smallest ratio 1.0136); splitting exact to 1e−14; A monotone → `indep_TA.json` | 2 s |
| `indep_S.py` | all 92 S leaves: the contour rebuilt from the stored design; s = ±(r − r_h) integrated as log s; q̂ at the 16 vertices of each leaf's parameter box. Closing boxes by the exact near-horizon map (no Taylor model), and the λ → 0 limit. W1 crossing heights | every certified bound ≤ the float minimum; closing boxes 0.4596 (+), 0.3970 (c) against certified 0.4147, 0.3534; crossings ≤ 0.0537 → `indep_S.json` | 3 s |
| `nh_check.py` | the closing box's near-horizon constants (m_f, the Rouché bound, M_G, M_W) by sampling; g_n by an FFT Cauchy integral; G(ρ₁) ≥ ζ₁² | all stored constants valid; g₁, g₂ agree → `nh_check.json` | 5 s |
| `indep_L.py` | Jost y at the stored L exits, by a Frobenius series in r (mpmath, 30 digits); the exit rule on 201 h² values per box. `indep_L.py out.json` checks 52 boxes at l ∈ {20, 21, 33, 100, 1000, 10⁵}; `indep_L.py out.json all` checks all 1760 boxes at l ∈ {20, 21, 100, 10⁵} | no violation in 1872 + 35 200 values; smallest margins 9.8e−4 (Re y₊), 5.6e−5 (Re y_c) → `indep_L.json`, `indep_L_all.json.gz` | 14 s; 266 s (all) |
| `bulk_and_affine.py` | the bulk facts (min W̃, U monotone, sympy identities); a randomized enclosure test of `p2.affine.AF` (it imports the class under test) | W̃ ≥ −0.0078910 > −99/12500; 160 020 checks, no violation → `bulk_and_affine.json` | 4 s |
| `indep_points.py` | W_l ≠ 0 at 150 random (l, ω), l ∈ [20, 2000], by a Riccati integration from both horizons; positive controls at the l = 20 and l = 50 PS modes | no zero; the secant converges to both modes (2.5e−8κ₋, 5e−7κ₋) → `indep_points.json` (and the 6-point trial `indep_points_test.json`) | 590 s |
| `k4_stops_indep.py` | D = \|y₊ − y_c\|/(\|y₊\| + \|y_c\|) at 13 of K4's 61 secant stops | D = 1 at all 13 → `k4_stops_indep.json` | 18 s |
| `cover_recheck.py` | the cover from the stored files with exact rationals: tilings, flags, overlaps, W1 and tail strips; 2000 random (l, ω) lookups up to l = 10⁶ | everything holds; nothing uncovered → `cover_recheck.json` | < 1 s |
| `pilot_compare.py` | for the top S leaves: the leading-order margin, the full q̂ and the certified bound (finding S5) | 0.025–0.032 → 0.015–0.018 → 0.009–0.014 → `pilot_compare.json` | < 1 s |

## Checks that run the author's code

| script | what it checks | result (at `0e329217`) | runtime |
|---|---|---|---|
| `replay_author.py` | exact replays: the three `+dirty` rows; the thin S leaves and both closing boxes; every T box and A in full | every produced field equals the certificate → `replay_author.json` | 20 s |
| `replay_all_SL.py` | `replay_all_SL.py S out.json`: all 92 S leaves. `replay_all_SL.py L out.json`: all 1760 L boxes | all identical → `replay_all_S.json`, `replay_all_L.json` | 205 s; 335 s |
| `closing_err_check.py` | finding S1: the closing box with the first version's real-segment error term and with a rigorous one (disc and complex-segment bounds) | margins unchanged to 8 digits; largest cell drop 8.2e−5 → `closing_err_check.json` | < 1 s |
| `k4_secant_trace.py` | the K4 float secant (`qnm_float`, dps 20) step by step at three stops (finding S2) | overshoot to Im ω ≈ −7…−13κ₋, return, zero step → `k4_secant_trace.json` | 194 s |
| `count_l20.py` | **backs the l = 20 clause of the P2-1 ledger row.** S3's reviewed counting pipeline (`p2.part3`, unchanged on the S4 branch) at l = 20 on R′₂₀ = [−a₂₀, a₂₀] × [−b, e], with a₂₀ from the certified P2-2 certificate `lg_PB_l20.json`. Usage: `count_l20.py OUTDIR [workers]` | winding 0 (ball 0 ± 7.4e−4), 223 segments, \|W₂₀\| ≥ 4.07 on the right half. With P2-2 and P2-0, W₂₀ has no zero in S_δ, independently of S4 → `count_l20/` (segment JSONL, summary, full pieces, launch record) | 2008 CPU-s (3 workers, ≈ 11 min) |

## The follow-up of 2026-10-19 (`followup_2026-10-19/`, at `21358225`)

| script / file | what it checks | result |
|---|---|---|
| `s1_recheck.py` | the fixed closing-box bounds (`nh_P_bounds`, `nh_P_disc`) against \|P_N\| in Arb at 4025 samples per side on the actual contour, and 642 tail samples; the first version's real-segment bound at the same samples; the new test's body with the old bound; both closing boxes recomputed | new bounds: no violation. Old lower bound: 813 violations per side (worst relative excess 2.8e−4 (+), 9.9e−6 (c)). The new test fails 46 times with the old bound, so it would have caught S1. Closing leaves recompute to the stored [0.41468628, 0.41468629] and [0.35342094, 0.35342095] |
| `compare_regenerated.py` | the regenerated certificates against the reviewed ones (read from git at `0e329217`) | S: 92 leaves, no value change; 273 printed intervals re-rounded outward, each by at most one unit in the 8th significant digit (≤ 1e−7 relative, ≤ 1e−5 absolute). L: leaves identical; summaries from ball endpoints. Controls: N4's strings re-rounded, timings |
| `replay_all_S.json`, `indep_S.json`, `indep_TA.json`, `cover_recheck.json`, `nh_check.json` | the earlier checks re-run against the regenerated certificates | the 92 S leaves replay identically (203 CPU-s); every independent check holds |
| `closing_err_check.json`, `pilot_compare.json` | re-runs (S1, S5) at the new head | unchanged, apart from the author's outward printing of the stored bounds |
| `pytest_fast.txt`, `pytest_slow.txt` | the whole P2 suite, fast and with P2_SLOW=1 | fast: 86 passed, 17 skipped (32 s); P2_SLOW=1: 103 passed, 1 pre-existing warning (660 s) |

## Notes

- **`indep_L_all.json.gz`** is `gzip -9 -n` of the 7 MB `indep_L_all.json`. The sha256 of the
  uncompressed file is `6f4eae30f15022ba87c893b0faec61715b0d3b58631399ba60e04c2fb12e4de3`.
- `count_l20/launch_*.json` records the scratch directory the run used. `count_l20.py` writes into the
  directory given on its command line and never touches `problems/P2/results` or `certificates`.
