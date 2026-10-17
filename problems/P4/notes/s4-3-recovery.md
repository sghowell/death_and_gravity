# S4-3 recovery: what survived the 2026-09-07 reboot, and how to finish

Status 2026-10-16 (recovery stage, branch `p4/theorem-c-s4-3`).  The
commutator strip exclusion is not certified.  `certificates/theorem_c/strip/tc_commutator.json`
was never assembled, every September region output was lost, and 1 of the 576
region pieces has a stored result (`parts/band/band-000.json`, the driver smoke
test).  Each of the following is backed by an artifact produced today from the
committed code:
* the exact identities E1/E2/E3/E5 (`certificates/theorem_c/checks.json`);
* Lemma D on a background that encloses the certified EC solution
  (`certificates/theorem_c/lemma_d_certified.json`).  The as-committed run
  (`parts/lemma_d_midpoint.json`) reproduces the old record exactly, but its
  centre series sits on the a3 midpoint (§3).
The four tube caches are rebuilt and pinned by content hash.  A persistent
driver, `scripts/s43_parts.py` (commit `76c07e8b`), replaces the lost scratchpad
drivers.  Paths below are relative to `problems/P4/` unless they start with
`problems/`.

**Update, 2026-10-16 evening (§11).**  The centre-ball defect of §3 is fixed in
the code (commit `ed4cab14`): `tc_commutator._centre_expansion` and
`tc_monotone.run_centre` now build the centre series on the certified (a*, μ*)
balls.  On the fixed code (commit `02e93bac`):
* P4-3 still holds (`certificates/theorem_c/tc_monotone.json` regenerated);
* Lemma D reproduces `certificates/theorem_c/lemma_d_certified.json` exactly;
* E1/E2/E3/E5 still hold.

The §8 (iii) centre variant fails (commit `ad16bbac`).  310 node-independent
region pieces and the §8 (ii) benchmark are running detached on the Mac; §11.6
gives PIDs, logs, expected finish and the stop command.  The commutator strip
exclusion is still not certified.

Sources for everything stated about September: the S4-3 agent's transcript
(`agent-a7accaae8fcc3cf60.jsonl`, 1671 records, 2026-09-01T01:12Z to
2026-09-05T18:40Z) and the coordinator's transcript (status through
2026-09-07T15:25Z), both under
`~/.claude/projects/-Users-seanhowell-dev-research-death-and-gravity/cf5d2a12-…`.
Transcript line numbers are given as `L<n>`.

## 1. What survived

* Commit `14da7d1e` (this branch): `src/p4/validated/tc_commutator.py` (1688
  lines), `tc_tm.py` (537), the Lemma D wiring in `tc_energy.py`,
  `tests/test_tc_commutator.py`, `notes/s4-commutator.md`, and
  `certificates/theorem_c/tc_energy.json` with the `lemma_d` block (ok, 61.6 s,
  nsub = 1100).  The committed `tc_commutator.py` is the final working-tree
  state: it contains the strip-depth fix (`if sdepth < 9`, line 1468; the
  09-04T01:22Z edit, L1529) and the progress logging added 09-05T01:13Z (L1621,
  L1624).
* The two transcripts, which hold every driver, launch command and logged result.
* The cloud node `ubuntu@216.86.163.117`, state unknown (Branch A, §6).

## 2. What was lost

The scratchpad `/private/tmp/claude-501/…/cf5d2a12-…/scratchpad/`:

* tube caches `tube_full.json` (344 steps), `s43/tube_gap.json` (107),
  `s43/tube_deep.json` (129), `s43/tube_band.json` (60);
* every region output: `centre_rep.json`, `tube_w0–w5`, `gap_w0–w5`,
  `gapd_w0–w3`, `sonicf_w0–w3`, `tubeband`, `tband2_a/b`, `tband3`,
  `tchunk_c0–c15`, `dchunk_d0–d7`, `centre_up_t1`, and all logs;
* the drivers `worker.py`, `gap_worker.py`, `gapd_worker.py`, `tubefix_worker.py`,
  `sonic_worker2.py`, `run_centre*.py`, `tchunk.py`, `dchunk.py`, `tband2/3.py`,
  `tubeband.py`, `build_gap_tube.py`, `build_deep_tube.py`, `curate_merge.py`,
  `dbg_deep1/2.py`, and the cloud's `/tmp/cworker.py`.

Nothing was merged: `curate_merge.py` (L1389) was written to run once the deep
band finished, and the deep band never finished locally.  The cloud production
round (16 `cworker` chunks, launched 2026-09-04T01:30Z, L1547) was never
retrieved; at the last check (2026-09-06T18:32Z) 2 of 16 chunks had finished,
both passing: `cw A1 ok= True worst= -63145.02 t=105341s`, and one more that
was not identified.

## 3. Re-verification today (committed code, Apple M4 Pro, flint 0.9.0)

| item | command | outcome | wall | artifact |
|---|---|---|---|---|
| E1, E2, E3, E5 | `scripts/s43_parts.py checks` | E1 True; E2 (pmin = 14, ok) True; E3 (True, True, True); E5 (True, True, True) | 48 s (E5 42 s) | `certificates/theorem_c/checks.json` |
| Lemma D as committed | `s43_parts.py lemma-d --centre-balls midpoint` | ok; 1100/1100 subranges; q2 ≥ 1.9618364, mB ≥ 3.4752564, mC ≥ 15.348300, pole gaps ≥ 0.99952, 0.99999; **identical to tc_energy.json** | 62.6 s | `parts/lemma_d_midpoint.json` |
| Lemma D, certified background | `s43_parts.py lemma-d --centre-balls certified` | ok; 1100/1100; q2 ≥ 1.9618364, mB ≥ 3.4752428, mC ≥ 15.348292, pole gaps ≥ 0.99952, 0.99999 | 61.5 s | `certificates/theorem_c/lemma_d_certified.json` |
| tube caches | `s43_parts.py build-caches --detach` | full 344 steps (187 s), gap 107 (121 s), deep 129 (8 s), band 60 (36 s): same step counts and step grids as September; content hashes reproducible (deep, band rebuilt twice: identical) | 3.2 min | `strip/caches/*.json` |
| test_tc_reduction.py | pytest | 9 passed | 23 s | — |
| test_tc_energy.py (+ `P4_TUBE_CACHE`) | pytest | 6 passed, 1 skipped (P4_SLOW); with P4_SLOW=1 the field identity passed | 12 s + 101 s | — |
| test_tc_commutator.py, fast set | pytest | 8 passed (incl. `test_sonic_closing_direct_phase`, 68 s), **2 failed** (`test_commutation_value_consistency`, `test_kdiv_roundtrip`), 1 skipped (no tc_commutator.json); **`test_centre_first_piece` not completed: stopped after 15 min (924 s wall)** | 104 s + 68 s + 924 s | — |
| test_tc_commutator.py, P4_SLOW | pytest | `test_e2_layer_identity` failed as written (14 ≥ 16); after the assertion fix (§9) it and `test_lemma_d_rerun` passed; `test_sonic_laurent_piece_full` passed | 71 s, 146 s | — |

`test_centre_first_piece` is one centre piece, y ∈ [1e-8, 1.6e-8], with 2 of
the 12 ε-cells.  Its docstring promises "~2–3 min" for the whole file, but the
test alone had not finished after 15 minutes and was stopped under the 15-minute
rule.  This is consistent with the centre cost in §4: no committed-code centre
run has ever finished.

Driver smoke test (one cheap piece, through the persisted driver):
`work --ids band/band-000` gave pass-both, margin −3.65e7, 77 s.  Then
`replay band/band-000` reproduced the margin bit for bit.  Artifacts:
`parts/band/band-000.json` and `parts/replays/band/band-000.json`.  This is the
first stored region piece; 575 remain.

The two fast failures are tightness assertions, not inconsistencies: in both,
the difference ball contains 0 (`|B1 − B1c| ⊂ [± 2.96e8]` against |B1| ≈ 1.2e9;
`|lhs − rhs| ⊂ [± 2.7e4]` against |lhs| ≈ 6.3e10), but the tests demand relative
widths below 1e-4 and 1e-8, which the interval hulls on ξ ∈ [8e-7, 8.2e-7] do not
reach.  The same two tests failed in every September run (`...s...FF.`,
pytest_final.log 09-02, pytest_final2.log 09-04; the agent loosened the first one
from 1e-6 to 1e-4 at L1051 and it still failed).  They are left unchanged; S4-4
should re-base them on point ξ so that the balls are tight.

Why Lemma D was run twice: `tc_commutator._centre_expansion` (and
`tc_monotone.run_centre`) build the centre series from the midpoint of
`results/a3_midpoint.json` only (radius ≤ 1e-41), without the stored Krawczyk
offset `K_minus_m`.  That point is not the EC solution: the certified balls are
a* = a_m + 7.031e-12 ± 8.87e-16 and μ* = μ_m + 2.3757e-9 ± 4.61e-14
(`modecount.certified_centre` and every tube build use these).  The committed
Lemma D record therefore certifies a background whose μ differs from the EC
value by 2.4e-9 (and a by 7e-12).
With the certified balls (`centre_balls="certified"`, the driver's default for
Lemma D and for the centre region) Lemma D passes with essentially unchanged
margins, so Lemma D now rests on a correct artifact.  The same correction is
needed in `tc_monotone.run_centre` (P4-3's centre tail) and should become the
code default in S4-4 (§9).  [Done later on 2026-10-16: both functions build
the certified balls in the committed code; see §11.1–§11.3.]

## 4. Inventory of the region parts

All pieces are listed in `certificates/theorem_c/strip/plan.json` (576 pieces, 33
excluded with reasons).  Each piece is computed by the committed region function
with a one-piece slice; the call is recorded in the part (`extra.region_call`).
Driver: `PY=.venv/bin/python; $PY problems/P4/scripts/s43_parts.py`.

| part | pieces | inputs | region call (per piece) | last September status | September cost (logged) | estimate, Mac P-core h |
|---|---|---|---|---|---|---|
| checks E1/E2/E3/E5 | 1 | — | `checks` | claimed in note §7, no artifact | — | done today, 0.01 |
| Lemma D | 1 | a3_midpoint.json | `lemma-d` → `run_lemma_d(at)` | recorded in tc_energy.json | 61.6 s | done today, 0.02 |
| centre `centre/c00–c29` | 30 | centre series (certified balls) | `centre_piece(at, ce, lo, hi, 14, eps_cells_default(), 320)` = `run_centre` body; y-pieces `y_subranges(1e-8, 1.6)` | `centre_rep.json` 30/30, worst −5.99, 1022 s — **produced 09-01T04:18Z by pre-direct-phase code** (τ-quadratic from τ = 14; `direct_S_boxes` was added at 04:32Z, L492); committed-code runs never finished: 06:23Z run killed after 26.0 h wall / 25.3 h CPU (L1364, L1372), 09:45Z run killed after 8.4 h (L1331); today a 2-cell test piece did not finish in 15 min | ≥ 25.3 CPU-h, unfinished | 45–90 (≥ 25 measured; the §8 (iii) minutes-long variant fails, §11.5) |
| tube idx 0–45 (`tube/full-000…045`, x ∈ [−0.267, −0.05]) | 42 | tube_full | `run_tube(at, tube_full, 14, x_min=−4.5, x_max=−0.05, step_slice=(i, i+1), degrade=True)` | tube_w0 (code of 09-01T06:58Z, no degrade): 42 pass, 4 fail with `TM inverse: \|q\| hull >= 1` at idx 29–32 (excluded, → band) | 28228 s wall / 46 steps | 4–7 |
| band `band/band-000…017` (x ∈ [−0.17, −0.134], h = 0.002) | 18 | tube_band | `run_tube(at, tube_band, 14, x_min=−0.1705, x_max=−0.134, …, degrade=True)` | tubeband 16/16 + tband3 2/2 pass (degrade runs; worst −2.6e7, −2.5e7); today `band/band-000` passed both phases | 1294 s incl. tube build | 0.5 |
| tube idx 46–91 (x ∈ [−0.784, −0.267]) | 46 | tube_full | as above | all pass, worst −1.14e6 | 8434 s | 1.5–2.5 |
| tube idx 92–137 (x ∈ [−1.349, −0.784]) | 46 | tube_full | as above | all pass, worst −1.69e5 | 3856 s | 0.6–1 |
| tube idx 138–166 (x ∈ [−1.524, −1.349]) | 29 | tube_full | as above | pass (inside tube_w3) | tube_w3: 59708 s / 46 steps incl. 12 failures | 1–5 |
| tube idx 167–183, "band A" (x ∈ [−1.814, −1.514]) | 17 | tube_full | as above | old cap (sdepth < 3): 12 fail in the direct phase at (τ ≈ 14, ε ≈ 0, first ℓ-strips); tband2_a/b (degrade) fail, +170.6 / +104.7, 34 h / 27.4 h; **cloud A1 (idx 167–174, sdepth < 9): pass, worst −63145, 105341 s on one vCPU**; A0 (idx 175–183): unknown | 3.7 vCPU-h/step (A1) | 30–65 |
| tube idx 184–244, "band B" (x ∈ [−3.05, −1.814]) | 61 | tube_full | as above | old cap: tube_w4/w5 (idx 184–273) killed after 25.4 h wall (24.7 CPU-h each), unfinished; tchunk c0 fail +209 (39.8 h), c2 +259.5 (57.7 h), c3 +824.9 (78.4 h), others killed 09-05; cloud B0–B7 (8, …, 5 steps): not finished at 65 h | ≥ 8 vCPU-h/step (cloud); 5–8 CPU-h/step (old cap, Mac) | 245–490 |
| deep `deep/deep-000…128` (x ∈ [−4.56, −2.95], h = 0.0125) | 129 | tube_deep | `run_tube(at, tube_deep, 14, x_min=−4.6, x_max=−2.9, …, degrade=True)` | dchunk d0–d7 (old cap) stopped/killed unfinished; cloud C0–C5 (21–22 steps): not finished at 65 h (at most one) | ≥ 3 vCPU-h/step | 195–520 |
| gap `gap/gap-000…106` (x ∈ [−0.055, −0.0085]) | 107 | tube_gap | `run_tube(at, tube_gap, 14, x_min=−0.056, x_max=−0.008, …, degrade=True)` | τ-quadratic phase fails on idx 0–88 (ξ ≤ 0.0357; bounds +8e8 … +7.8e10); gapd_w0–w3: idx 0–89 pass the direct phase (worst −0.072 … −0.78), τ > 224 open; gap_w5 idx 90–106 and idx 89 pass both phases (worst −2.96e7) | 92611 s wall | 15–25 |
| sonic `sonic/s00–s50` (ξ ∈ [0, 0.0091]) | 51 | sonic series | `run_sonic(at, 14, piece_slice=(i, i+1))` | sonicf_w0–w3 (code of 09:45Z): 51/51 pass the direct phase (worst −5.9e−5), τ > 224 open | 16062 s wall | 3–4 |
| excluded: tube idx 245–273 (h = 0.05, x ≤ −3.05) | 29 | — | — | tchunk c11–c15: `TM inverse` failure in ~12 min each → replaced by deep | — | — |
| excluded: tube idx 29–32 | 4 | — | — | `TM inverse` at h = 0.0075 → replaced by band | — | — |

Total ≈ 540–1200 Mac P-core hours, of which the deep band x ∈ [−4.56, −1.514]
(band A + band B + deep, 207 pieces) is 470–1070 and the centre 45–90.

How the estimates were made.  September wall times on the Mac were taken under
1.5–2× oversubscription (15–34 concurrent workers on 12 cores; the logged CPU
shares were 60 % at 09-01T09:57Z and 94–100 % later), so wall × (0.6–0.9) is used
as CPU time.  Cloud times are per dedicated vCPU: the coordinator's "32 procs on
16 vCPUs, 2× oversubscribed" counted the idle `uv run` parent of each worker;
`uptime` gave load 16.00 (L1579) and `ps` showed each python worker with CPU time
equal to elapsed time (L1606: `01:42:30` CPU in `01:42:31`).  A cloud vCPU is
taken as 0.5–1 Mac P-core (the coordinator's "2× slower" was never measured).
The centre row is a lower bound only: no committed-code centre run ever finished.
The code version matters for every September number except the cloud round:
tube_w0–w5 ran the code of 06:58Z (before the 08:33Z branch-R sign fix and the
08:42Z offset-aware direct chains), centre_rep the code of 04:18Z, gap 09:23Z,
sonicf 09:45Z, gapd/tubeband/tband after 18:12Z (with `degrade`), and only the
cloud round (and nothing local) ran with the strip-depth fix.

Corrections to the note found in the logs: the near-sonic open band was
ξ ≤ 0.0357 (gap idx ≤ 88) in the September data, not "~0.048"; the main tube is
241 + 18 (band) + 129 (deep) pieces, not 274 steps of tube_full; the gap tube
passed both phases on 18 of its 107 steps, not "107/107 both phases".

## 5. The driver (`scripts/s43_parts.py`) and the result layout

* `build-caches [--detach]`: the four September tube constructions, verbatim
  (recipes recorded in `strip/caches/<name>.json`); caches go to
  `problems/P4/data/cache/` (gitignored by the root `.gitignore`), pinned by a
  content hash over (K, steps) that excludes the wall-clock build time.  `work`
  refuses a cache whose content hash differs from its record.
* `plan`: enumerates every piece (`strip/plan.json`, committed).
* `work` / `launch --workers N --nice 10`: per-piece results in
  `parts/<region>/<id>.json`; a worker claims a piece with an `O_EXCL` lock file
  (stale locks of dead local workers are reclaimed), writes atomically, and skips
  pieces that already have a valid result (spec equal to the plan entry, cache
  content hash equal, every recorded source file still at its recorded sha256;
  certified failures are also skipped unless `--retry-failed`).  `launch` starts
  detached workers (new session, SIGHUP ignored, `nohup`, `nice`) with logs and a
  launch record in `problems/P4/logs/s43/` (gitignored).  All paths derive from
  the script's location, so the same command runs on the node.
* Selection for machine splits: `--regions`, `--ids`, `--ids-file`,
  `--exclude-ids-file`, `--x-range LO,HI`, `--shard K/M`; `list` prints the ids.
* `status`, `replay ID` (recompute one piece, compare, never overwrite),
  `import-cloud` and `cloud-manifest` (Branch A), `merge`.
* Added 2026-10-16 evening (commit `ed4cab14`).  Plan entries carry `tq`
  (False = direct phase only; gap idx 0–88 and all sonic pieces, set by
  `REGIONS[...]["direct_only"]`) and, when set, `tstar_fac`.  Both are part of
  the spec, so a direct-only part never stands in for a both-phase one.
  New commands:
  * `lemma-d --replay` runs the committed `run_lemma_d` and compares the result
    with the stored record;
  * `p43-centre` runs P4-3's centre tail and replays the old midpoint record;
  * `bench` times `analyze_piece` with `tq=True` and `tq=False` and records the
    phase split.  It is detachable and resumable; results go to
    `results/theorem_c/bench/`, which `merge` does not read.
* `merge` checks coverage: the union of the x-intervals of passing pieces must
  contain (−∞, 0] (centre: y-pieces mapped to x = ln y; sonic: x = −ξ); pieces
  with `open` count only for τ ∈ [14, 224].  It writes
  `certificates/theorem_c/strip/tc_commutator.json` only when every planned piece is
  present and passing, the checks are valid for the current code, and no
  unaccepted code variant is mixed in; otherwise it writes
  `tc_commutator_partial.json` with the missing and failed lists.

Why the per-piece parts are committed: they are the certificate's evidence, each
is about 4 KB (≈ 2.3 MB for all 576), and committing them as they land is what
would have saved September's work.  The merged certificate is regenerated from
them deterministically; logs, locks and tube caches stay out of git.

## 6. Branch A: the node still has the September results

All `ssh`/`scp` below are for Sean to run; this stage did not contact the node.

1. Inventory:
   `ssh ubuntu@216.86.163.117 'uptime; ls -la ~/out ~/caches; grep -h "^cw " ~/out/cw_*.log | sort'`
   Expect 16 lines `cw <tag> ok= … worst= … t=…s` (A0, A1, B0–B7, C0–C5) and 16
   `cw_*.json`.  Missing tags mean their chunks died; those pieces go to Branch B.
2. Code identity: run
   `ssh ubuntu@216.86.163.117 'cd ~/death_and_gravity && sha256sum problems/P4/src/p4/__init__.py problems/P4/src/p4/validated/*.py' > node.sha`
   and compare with the local tree (`sha256sum` of the same files in this
   worktree).  Expected: identical.  The node's `tc_commutator.py` on disk was
   overwritten at 09-05T01:14Z with the committed version, but the cworker
   processes (started 09-04T01:30Z) executed the version scp'd at 09-04T01:22Z,
   which is the committed file minus the two `P4_PROGRESS` print blocks.  Its
   sha256 is `0c5a8e67ce9f2661475f8833fbbc71957c9d9a1238f9ca02898accf16f042c7f`;
   `s43_parts.py cloud-manifest --out …` regenerates it and writes the variant to
   `problems/P4/logs/s43/` for `diff`.  The committed manifest is
   `results/theorem_c/parts/imported/cloud-2026-09/code_manifest.json`.
3. Cache identity: on the node,
   `cd ~/caches && python3 -c "import json,hashlib; [print(f, hashlib.sha256(json.dumps(dict(K=d['K'], steps=d['steps']), sort_keys=True, separators=(',',':')).encode()).hexdigest()) for f in ('tube_full.json','tube_deep.json') for d in [json.load(open(f))]]"`
   Expected: `tube_full.json cb2531921be5ff5c387f2f1b424f351c7fa378f80f71c97dc7edd464de0b1b6e`,
   `tube_deep.json ef544d927ef12534ed38a4f04044a4bbb6d61f5a45f0f7c5e20ec00775cfa493`
   (today's rebuilds).  If they differ, the node's steps are not the planned ones:
   copy the node caches here, re-plan on them, and import against that plan.
4. Pull back (raw reports are provenance, commit them):
   `mkdir -p problems/P4/results/theorem_c/parts/imported/cloud-2026-09/raw && scp 'ubuntu@216.86.163.117:out/cw_*' problems/P4/results/theorem_c/parts/imported/cloud-2026-09/raw/`
5. Import: `$PY problems/P4/scripts/s43_parts.py import-cloud problems/P4/results/theorem_c/parts/imported/cloud-2026-09/raw --code-manifest problems/P4/results/theorem_c/parts/imported/cloud-2026-09/code_manifest.json --tube-hash full=<hash from 3> --tube-hash deep=<hash from 3>`.
   Mapping: A0 → `tube/full-175…183`, A1 → `tube/full-167…174`, B0–B7 →
   `tube/full-184…244`, C0–C5 → `deep/deep-000…128`, matched by (cache, x, h);
   unmatched steps are reported.
6. Spot replays (fresh compute, ≈ 2–4 h each on one core): `replay tube/full-167`
   and `replay deep/deep-064`; identical margins are expected because the code is
   math-identical and the tubes are identical.
7. Remaining pieces: centre (30), tube idx 0–166 (163), band (17 left), gap
   (107), sonic (51): 368 pieces, ≈ 70–135 core-h, of which the centre is
   two thirds.
   On the Mac: `$PY problems/P4/scripts/s43_parts.py launch --workers 8 --nice 10 --accept-code tc_commutator.py=0c5a8e67ce9f2661475f8833fbbc71957c9d9a1238f9ca02898accf16f042c7f`
   (imported parts are then recognised as valid and skipped).  Wall ≈ 9–17 h,
   set by the slowest centre pieces (the §8 (iii) centre variant fails, §11.5).
   [Update 2026-10-16 evening: tube idx 0–137, band, gap and sonic are running
   (§11.6); afterwards the node-independent remainder is the centre (30) and
   tube idx 138–166 (29, in no cloud chunk: A1 starts at idx 167), held back
   with the rest of idx 138–244 until the node check.]
8. Merge: `$PY problems/P4/scripts/s43_parts.py merge --accept-code tc_commutator.py=0c5a8e67… --accept-reason "cloud chunks 2026-09-04: executed tc_commutator.py = 14da7d1e minus the two P4_PROGRESS blocks"`.
   The τ > 224 band stays open in `scope` until S4-4.

## 7. Branch B: the node is gone or its results are incomplete

0. Prerequisite experiments (≈ 1 day, decide before the big launch): the centre
   variant and the τ-quadratic profiling of §8.
1. Caches on each machine: `build-caches` (≈ 3 min); `work` refuses any cache
   whose content hash differs from the committed record.  On the node, check out
   this branch first (its repo copy predates the S4-3 files).
2. Split (deep band shared, everything cheap on the Mac):
   `list --x-range=-4.6,-1.52 --shard 0/2 > node.ids` (104 pieces);
   node: `launch --workers 16 --nice 10 --ids-file node.ids`;
   Mac: `launch --workers 8 --nice 10 --exclude-ids-file node.ids` (472 pieces;
   workers take pieces in order of estimated cost, so the cheap regions finish
   first and the Mac's 103 deep-band pieces come last).  Use 10 Mac workers if
   the machine is otherwise idle (4 of the
   12 cores are efficiency cores, worth ≈ 0.3–0.5 of a P-core each).
   Pull the node's `parts/` back with `rsync` at least daily and commit.
3. Wall clock: total ≈ 540–1200 P-core h against ≈ 16 P-core equivalents (8 Mac
   P-cores + 16 vCPUs at 0.5) gives 34–75 h with perfect balance.  With the 1:1
   split of the deep band, the node finishes its 235–535 h half in 30–67 h if a
   vCPU is half a Mac core (15–34 h if it is a full one).  The Mac finishes its
   70–135 h of cheap pieces plus 235–535 h of deep band in 38–84 h.  Plan on 2–4
   days; single deep-band pieces of up to ≈ 8 core-h set the tail.  The Mac
   alone (8 workers) needs 68–150 h.
4. Merge as in Branch A step 8 (no `--accept-code`).

## 8. The large-|κ| option and the S4-4 route

The option (note §8 (ii), §10): exclude τ > 224 by a large-|κ| perturbative
argument, making both the τ-quadratic phase and the near-sonic sliver
unnecessary.

**Feasibility.**
1. Whatever replaces the τ-quadratic phase must hold uniformly up to the sonic
   point, and that is exactly where large-τ expansions stop being uniform.  On the
   band the commuted condition scales as S ℓ²/τ² = F(ℓτ, ε) + o(1) (September
   float check, `corner3.py`: F(0+) = −0.90, F(1) = −3.21, F(∞) = −15.5 at ε = 0;
   F ≤ −41 at ε = 3; ℓ-independent to three digits between ℓ = 1e-4 and 1e-5).
   The natural variable is the sonic-layer variable s = ℓτ, i.e. the layer
   |x| ≲ 1/|κ| of theorem-c-plan.md §2.3.  A perturbative exclusion that treats κ
   as the only large parameter must resolve that layer; it does not go around the
   obstruction that defeats the τ-quadratic phase.
2. In the matching-function form (route (a) of the plan, ledger P4-7) the input is
   the law E(κ) = C κ^{−1} e^{Φ₊κ}(1 + d/κ + O(κ^{−2})), C ≈ 0.633, d ≈ −8.5.  It is
   CONJECTURED: float fits on three rays, |κ| ∈ [30, 300].  Making it an exclusion
   needs Lemma L1 (the sonic-layer connection coefficient c_s, uniform in
   arg κ ∈ [−π/2, π/2]), which the plan rates "genuine research", plus Levinson
   bounds with explicit constants (L2) and the centre layer (L3).  To retire the
   τ-quadratic phase the certified threshold must satisfy K₀ ≤ 224.  The error
   factor has the form (C₁/|κ|)e^{C₂}, with C₂ a weighted integral of the coupling
   matrix B̃ over the Levinson zone.  Integrals of B̃ there are of size 1–10: the
   diagonal one that enters the leading factor is already J = 6.27 (e^J ≈ 530).
   So a rigorous K₀ ≤ 224 needs near-optimal constants.  If K₀ comes out above 224, the
   direct phase must be extended to K₀, which costs more compute, not less; the
   route-(c) contour bridge is comfortable only up to |κ| ≈ 130 (plan §5 (v)).
3. In the energy framework the same layer is a bounded computation: the joint
   (w = 1/κ, ℓ) Taylor-model engine of note §8 or, more simply, a scale-invariant
   layer bound S ℓ²/τ² ≤ F(ℓτ, ε) + error, with F built from sonic hats that the
   code already certifies (u2, u1, tt, w0; b0^(−3) = 2.10/ξ³).  The float margin is
   at least 0.9, the work was estimated at about 150 lines on top of `tc_tm.py`,
   and it touches the 51 sonic and 89 gap pieces only.  Its w-pseudo-box
   (w ∈ [0, 1/224]) applied to the direct chains away from the layer could also
   replace `sup_ratio_boxes`, retiring the τ-quadratic phase without route (a).
4. Cost evidence from September.  On the gap tube the τ-quadratic remainder
   machinery was roughly two thirds of a full first attempt.  gap_w0–w4 (attempt
   only) took 61910 s for 90 steps; gapd (attempt plus direct-only rerun) took
   81089 s.  So the direct phase alone was about 19179 s.  This is a wall-time
   difference between runs under similar but not identical load.  On those steps
   the τ-quadratic phase then failed, so that compute was wasted.  On the centre the opposite holds: the
   τ-quadratic-only code of 04:18Z certified all 30 pieces in 1022 s, while
   the committed direct phase on τ ∈ [14, 224] has not finished in 25 CPU-h.  On
   the deep band (90 % of the cost) the split is unknown.

**Recommendation.**  Do not make route (a) the S4-4 route, and do not hold the
re-run for it: P4-7 is a conjecture, L1 is open research, and an explicit
K₀ ≤ 224 is unlikely.  Do settle the τ > 224 strategy before spending compute on
the τ-quadratic phase, because the direct phase on τ ∈ [14, 224] is needed under
every route and the τ-quadratic phase is not.  Concretely, S4-4 should:
(i) build the energy-framework layer engine for the band ξ ≤ 0.036, where the
τ-quadratic phase fails anyway, and run those 140 pieces direct-only (a `tq=False`
driver mode, about 4× cheaper than `degrade` there);
(ii) before launching the deep band, time `analyze_piece` with `tq=True` and
`tq=False` on two deep steps (`tube/full-200`, `deep/deep-064`; about one day on
2 cores).  If the τ-quadratic phase is a minor share, keep the committed
both-phase runs, which passed outside the band in September.  If it dominates,
run direct-only everywhere and let the layer engine's w-pseudo-box carry τ > 224;
(iii) try the centre with the τ-quadratic phase from τ = 14 (`TSTAR_FAC = 1`
for the centre pieces, recorded in each part).  One piece takes minutes, and if
`sup_ratio_boxes` closes from κ = 14 as it did in the 04:18Z code, the centre
drops from ≥ 25 core-h to minutes.
Route (a) stays the long-term fallback and a useful float diagnostic.

**Outcomes, 2026-10-16 evening.**
* (iii) was tried and fails, quickly (§11.5).
* (i): the direct-only pieces are running (§11.6).
* (ii): the benchmark is running.  Its logs already show the τ-quadratic
  bounds taking about 4 s per piece; the direct phase was still running at
  hand-off (§11.6).
* Update 2026-10-17 (`s4-4-layer-engine.md` §§9–13): the S4-4 layer engine,
  with a per-piece margin scale, polynomial-level derivative blocks and a
  centred enclosure (`tc_layer_mv`), certifies tube/full-170, -200, -240 and
  deep/deep-000, -064, -128 for τ ≥ 224 and for τ ≥ 14 (both phases at once)
  at half of each piece's float scale, in 1–4.6 min per piece.  Projected
  cost of tube 138–244 plus deep 0–128 that way: about 8 core-h instead of
  470–1070.  Using it needs a reviewed merge rule (a passing `layer_t14` part
  in place of a missing S4-3 part); not implemented.

## 9. Other findings to carry forward

* `tc_commutator._centre_expansion` and `tc_monotone.run_centre` use the
  a3 midpoint without `K_minus_m` (§3).  Fix: build (a, μ) as
  `s43_parts.certified_centre_balls()` does, then re-run Lemma D (62 s, done today
  in the driver) and the P4-3 centre tail (`tc_monotone`, seconds).
  [Done 2026-10-16 evening, §11.1–§11.3.]
* `run_lemma_d` raises Z to offset 16 trusting E2, but does not check that the
  returned pmin satisfies pmin + 2 ≥ 16 (it is 14, so the run is sound); add the
  assertion.  [Done 2026-10-16 evening: it raises unless pmin + 2 ≥ 16, §11.1.]
* `test_e2_layer_identity` asserted `val >= 16` on pmin (= 14) and failed; fixed
  to `val + 2 >= 16`, which matches `check_layer_identity`'s contract and the
  note.  It now passes.
* `run_tube(degrade=True)` recomputes the whole direct phase after a failed
  τ-quadratic attempt; on the band that doubles the direct cost.  In the
  `tq=True` call, `analyze_piece` has already run `direct_S_boxes` on every
  ε-cell with the same arguments as the `tq=False` call.  S4-4 could
  therefore reuse those results.  A rerun is needed only when the
  τ-quadratic setup raises before the cell loop (e.g. thU ≥ 1).
* The tube caches are ≈ 9 MB in all and rebuild in 3 min; they stay out of git and
  are pinned by content hash.
* Centre τ-quadratic structure (§11.5): on c29 the τ² coefficient K2 carries
  y² while the τ-linear remainder slack has a y-independent part (0.441).  A
  τ-quadratic phase that is to start near τ = 14 in the centre needs
  remainder bounds whose τ-linear part also vanishes like y² (or a
  y-dependent tstar).  The S4-4 layer engine should check this.

## 10. Draft ledger text (for Sean; CLAIMS.md not edited)

* P4-8: "conditional on Lemma D" can be dropped on the strength of
  `certificates/theorem_c/lemma_d_certified.json` (Lemma D with the certified
  (a*, μ*) balls, ok, 1100/1100 subranges, q2 ≥ 1.9618, mB ≥ 3.4752, mC ≥ 15.348;
  replay `s43_parts.py lemma-d`).  The midpoint-ball record in tc_energy.json
  reproduces exactly but does not enclose the EC background.
* P4-9 (commutator strip): no certificate yet; status after the re-run per §6/§7.

Revised drafts, 2026-10-16 evening (supersede the two items above):

* **P4-3**: the status stays CERTIFIED.  Suggested addition to the claim
  text: "(centre tail on the certified (a*, μ*) balls since 2026-10-16)".
  * Evidence: `certificates/theorem_c/tc_monotone.json`, regenerated on the fixed
    code.  The centre tail t ∈ (0, e^{−8}] is now evaluated on balls that
    enclose the certified (a*, μ*): a3 midpoint + Krawczyk offset,
    `p4.validated.ec_centre`.  The previous centre block used the bare
    midpoint, which does not enclose them.
  * All six t-factored margins stay positive:
    * v_rel/t ≥ 0.40433;
    * (v_rel)'/t ≥ 0.40433;
    * −ρ̂′/t² ≥ 82.744;
    * 1/3 − v_rel² ≥ 0.33333;
    * −t²D ≥ 1.5292;
    * t²D′ ≥ 3.0584.
  * None moved by more than 1.2e-7 from the midpoint run.  The sonic and tube
    (344 steps) blocks are unchanged bit for bit.
  * Provenance: `certificates/theorem_c/p43_centre.json`.  It also replays the
    old midpoint record exactly.
  * Replay: `cd problems/P4 && PYTHONPATH=src python -m p4.validated.tc_monotone
    data/cache/tube_full.json`, or `scripts/s43_parts.py p43-centre`.
* **P4-8**: suggested status "CERTIFIED (unconditional since 2026-10-16:
  Lemma D closed in S4-3)".
  * Lemma D (deep tail x ≤ −4.5) is certified with centre balls that enclose
    the certified (a*, μ*): `certificates/theorem_c/lemma_d_certified.json`.
  * It is reproduced exactly through the fixed committed code path:
    `results/theorem_c/parts/replays/lemma_d_certified.json`.
    `tc_commutator.run_lemma_d` there uses `ec_centre` and now also checks
    the E2 condition pmin + 2 ≥ 16.
  * Result: ok, 1100/1100 subranges; q2 ≥ 1.96184, mB ≥ 3.47524,
    mC ≥ 15.3483; pole gaps ≥ 0.99952, 0.99999.
  * Replay: `scripts/s43_parts.py lemma-d --replay`.
  * The `lemma_d` block in `tc_energy.json` is the old midpoint run.  It
    reproduces (`parts/lemma_d_midpoint.json`) but does not enclose the EC
    background; cite `certificates/theorem_c/lemma_d_certified.json` instead.
* **P4-9** (commutator strip): no certificate yet; status after the re-run
  per §6/§7/§11.

## 11. 2026-10-16 evening: integrity fix, re-runs, centre experiment, launch

### 11.1 Integrity fix (commit `ed4cab14`)

* New module `src/p4/validated/ec_centre.py` with two functions:
  * `certified_centre_balls()` returns the a3_midpoint.json midpoint plus its
    Krawczyk offset `K_minus_m`.  This is the driver's construction exactly.
    At 256, 320 and 384 bits the parsed balls contain the printed decimal
    balls (endpoint check, also done at run time, widening otherwise).
  * `certified_centre_expansion(prec, K)` builds the centre series on those
    balls.
* `tc_commutator._centre_expansion` (Lemma D, centre region) and
  `tc_monotone.run_centre` (P4-3 centre tail) now use the new module.
  * Test `test_centre_expansion_encloses_certified_background`: the series'
    ŵ = μe^{2a} and n̂ = e^{−a} contain those of modecount's certified balls.
    Its ŵ is disjoint from the midpoint value.
* `run_lemma_d` now raises unless the E2 layer minimum satisfies
  pmin + 2 ≥ 16 (test `test_lemma_d_requires_layer_valuation_16`).
* New options, defaults unchanged:
  * `run_tube(tq=False)` and `run_sonic(tq=False)` run the direct phase only.
  * `analyze_piece`, `centre_piece` and `run_centre` take `tstar_fac`.  With
    tstar_fac ≤ 1 the τ-quadratic phase covers all τ ≥ b1, and
    `direct_S_boxes` is not called because its range [b1, tstar] is empty.
  * `tq=False` with tstar_fac ≤ 1 is refused (test
    `test_direct_only_needs_a_direct_range`).
* Driver changes: §5.

### 11.2 Re-runs on the fixed code (commit `02e93bac`)

| item | command | outcome | wall | artifact |
|---|---|---|---|---|
| P4-3, all three regions | `cd problems/P4 && PYTHONPATH=src python -m p4.validated.tc_monotone data/cache/tube_full.json` | PASSES.  Centre tail on the certified balls: vrel_pos ≥ 0.40433444061, dvrel ≥ 0.40433436552, neg_drho ≥ 82.743524313, wG_pos ≥ 0.33333331494, neg_Ds ≥ 1.5291792671, dD_pos ≥ 3.0583589280.  Sonic minima and tube (344 steps) minima identical to the stored record; only `tube_info.time` (cache build time) differs | 2.4 s | `certificates/theorem_c/tc_monotone.json` (regenerated) |
| P4-3 centre tail, with provenance | `s43_parts.py p43-centre` | as above; the midpoint replay reproduces the old centre block exactly.  Certified − midpoint: neg_drho +1.19e-7, dD_pos −4.3e-11, neg_Ds −2.2e-11, vrel_pos and dvrel +2.8e-12, wG_pos 0 | < 1 s | `certificates/theorem_c/p43_centre.json` |
| Lemma D, fixed code path | `s43_parts.py lemma-d --replay` | identical to `certificates/theorem_c/lemma_d_certified.json` in ok, c0, ycut, layer identity (14, True), nsub, failures (0 of 1100) and all five minima: q2 ≥ 1.961836442292799, mB ≥ 3.4752428204250765, mC ≥ 15.348292102988353, poleA ≥ 0.9995211043091139, poleFN ≥ 0.9999925864848808 | 64.7 s | `parts/replays/lemma_d_certified.json` |
| E1/E2/E3/E5 (tc_commutator.py hash changed) | `s43_parts.py checks` | E1 True; E2 (pmin 14, ok); E3 and E5 all True | ≈ 48 s | `certificates/theorem_c/checks.json` (overwritten; valid for the current code) |

### 11.3 P4-3 verdict

CLAIMS row P4-3 still holds.  The defect affected only the centre tail
t ∈ (0, e^{−8}].  On the certified balls every t-factored margin stays
positive:
* the smallest, 1/3 − v_rel² ≥ 0.33333331, is unchanged;
* no margin moves by more than 1.2e-7.

The sonic and tube blocks never used these balls (the tube's centre extension
was built on the certified balls), and they reproduce bit for bit.  Draft
ledger text: §10.

### 11.4 Tests

* `test_tc_commutator.py`, run with P4_SLOW=1 and `test_centre_first_piece`
  deselected (more than 15 min on the default path; not run): 14 passed,
  2 failed, 1 skipped (no merged certificate), 314 s.
  * The 2 failures are the known tightness assertions
    `test_commutation_value_consistency` and `test_kdiv_roundtrip`: the balls
    overlap; unchanged since September (§3).
  * The passes include the three new tests, `test_e2_layer_identity`,
    `test_lemma_d_rerun` (63 s, fixed code path) and
    `test_sonic_laurent_piece_full` (138 s, `run_sonic` default path).
* `test_tc_energy.py`, run with `P4_TUBE_CACHE`: 6 passed, 1 skipped (the
  P4_SLOW field identity; `tc_energy.py` is unchanged), 12 s.  The passes
  include `test_monotone_sonic_and_centre` (the fixed `run_centre`).

### 11.5 Centre experiment, §8 (iii) (commit `ad16bbac`)

Setup: `centre/c29` (y ∈ [0.00694, 0.01111], the outermost piece) with
`tstar_fac = 1`, run through the driver on the committed code.

Result: certified failure after 39 s wall (39 s CPU).
* The τ-quadratic bound is positive on 10 of the 12 ε-cells: +77.05 on
  ε ∈ [0, 0.03], falling to +1.22 on [2.4, 3.4].
* Only ε ∈ [3.4, 6] closes (−1.06, −2.36).

Diagnosis (`scripts/s43_centre_tq_diag.py`, output stored):
* K2, the τ² coefficient, has y-offset 2.  Its y-scaled hull is −5.9e-4 …
  −1.6e-2 over the cells.
* The τ-linear remainder slack R1 has an offset-0 part, 0.441.
* K0 + R0 ≈ −5.36.
* So the vertex of K2τ² + R1τ lies at τ ≈ 14 … 374, and its supremum (up to
  82) exceeds |K0 + R0|.

Consequences:
* K2 carries y² and the slack does not, so pieces closer to y = 0 should do
  worse.  This is expected, not run.
* The committed τ-quadratic phase does not reproduce the 04:18Z September
  centre pass.

As instructed, no further centre pieces were tried and no committed-path
centre run was started.  The plan keeps the centre on the committed default
(TSTAR_FAC = 16, direct boxes on [14, 224]).  Artifacts, in
`parts/experiments/centre-tstar_fac-1/`: `c29.json` (the part, with its
timing block) and `c29_tau_quad_diag.json`.

### 11.6 What is running

Launched 2026-10-17T03:38Z (2026-10-16 20:38 PDT) on the Mac, at nice 10.
Each process is detached: its own session, parent PID 1, SIGHUP ignored, nohup.
Six processes in all.

**Region workers** (tag `s43-20261016`): PIDs 44147, 44148, 44149, 44150.
* Logs: `problems/P4/logs/s43/s43-20261016-w0.log` … `-w3.log`.
* Launch record: `problems/P4/logs/s43/launch-s43-20261016.json`.
* Selection: `--ids 'tube/full-0??,tube/full-1[0-2]?,tube/full-13[0-7],band/*,gap/*,sonic/*'`,
  310 pieces:

  | pieces | count | phases |
  |---|---|---|
  | tube idx 0–137 without 29–32 | 134 | both phases, degrade |
  | band 0–17 | 18 | both phases, degrade |
  | gap 89–106 | 18 | both phases, degrade |
  | gap 0–88 | 89 | direct only, `tq = False` |
  | sonic 0–50 | 51 | direct only, `tq = False` |

* Order: estimated cost.  Band and tube come first; the both-phase gap
  pieces come last.
* Results go to `problems/P4/results/theorem_c/parts/<region>/<id>.json`.
* Expected finish: the estimate is 21 core-h.  Band pieces took 103–115 s
  against an estimate of 108 s.  That gives 5–10 h wall, i.e. about
  2026-10-17 09:00–14:00Z (02:00–07:00 PDT).
* At hand-off (2026-10-17T03:49Z):
  * the band region is complete: band-000 … band-017 are all pass-both,
    margins −3.65e7 … −2.50e7, 103–115 s wall each, committed;
    `tube/full-093` and `full-094` are pass-both and tube pieces are running.
  * The validity check passed on the first parts: spec equals the plan entry,
    source hashes are current, git head is `ad16bbac`, src is clean.
  * band-000 recomputed on the fixed code reproduces the morning's margin
    −36528931.26271176 bit for bit.

**Benchmark** (§8 (ii); tag `bench-20261016`): PIDs 44319 (`tube/full-200`)
and 44321 (`deep/deep-064`).
* Each process runs `analyze_piece(tq=True)` and then
  `analyze_piece(tq=False)` on its piece: `run_tube`'s per-step body at
  prec 320, no degrade.
* Wall and CPU are recorded, plus the phase split: `sup_ratio_boxes` +
  `tau_quad_check` against `direct_S_boxes`, per call.
* Logs: `problems/P4/logs/s43/bench-20261016-tube__full-200.log` and
  `…-deep__deep-064.log`.  Launch record: `launch-bench-20261016.json`.
* Results: `problems/P4/results/theorem_c/bench/tube__full-200__tq.json`,
  `tube__full-200__direct.json`, `deep__deep-064__tq.json`,
  `deep__deep-064__direct.json`.  These are not parts; `merge` ignores them.
* Expected: about a day, finishing 2026-10-17 late to 2026-10-18.  No Mac
  timing of these pieces exists; September gives ≥ 8 vCPU-h per band-B step
  and ≥ 3 vCPU-h per deep step.
* Early log observation (not a result; the JSON will hold the numbers).  On
  both pieces, after an 83–85 s setup:
  * `sup_ratio_boxes` (U, then the six families) took 2 + 2 s;
  * `tau_quad_check` took under 1 s;
  * then the direct phase started.  It was still running at hand-off.

**Not launched**, per instruction:
* tube idx 138–244, which includes band A (167–183) and band B (184–244);
* deep 0–128.

These wait for the cloud-node check (§6).  Tube idx 138–166 is in no cloud
chunk (A1 starts at idx 167), so it is needed under either branch.  The centre
(30 pieces, committed path only, ≥ 25 core-h) was not started either.

### 11.7 Checking, stopping, resuming, committing

`PY=.venv/bin/python` (from the repository root).

* Status: `$PY problems/P4/scripts/s43_parts.py status`.  It prints counts per
  region, the finished bench records and the running bench modes; add `-v`
  for each piece.  Also:
  * `tail -f problems/P4/logs/s43/s43-20261016-w*.log problems/P4/logs/s43/bench-20261016-*.log`
  * `ps -o pid,etime,time,command -p 44147,44148,44149,44150,44319,44321`
* Stop all: `kill 44147 44148 44149 44150 44319 44321`.
  * A piece in progress is lost and is recomputed on the next launch.
  * Lock files of dead processes are reclaimed automatically.
* Resume: valid results are skipped.
  * Regions: `$PY problems/P4/scripts/s43_parts.py launch --workers 4 --nice 10 --tag s43-<date> --ids 'tube/full-0??,tube/full-1[0-2]?,tube/full-13[0-7],band/*,gap/*,sonic/*'`
  * Bench: `$PY problems/P4/scripts/s43_parts.py bench --ids tube/full-200,deep/deep-064 --modes tq,direct --detach --nice 10`.
    Finished modes are skipped.
* Do not edit `problems/P4/src` while workers run.  A part records the source
  hashes read from disk when it is written; a changed file makes the part
  stale.
* Commit parts as they land:
  `git add problems/P4/results/theorem_c/parts problems/P4/results/theorem_c/bench`.
  `band/band-000.json` from commit `fef526e6` was computed on the pre-fix
  code; this launch has recomputed it.
* After the run, `merge` writes `tc_commutator_partial.json`: the centre,
  tube 138–244 and deep are still missing, and τ > 224 stays open on the
  direct-only pieces.

---

*Paths (2026-10-17): the files behind CERTIFIED rows moved unchanged from `results/` to `certificates/` (`problems/P4/certificates/README.md`). The paths in this note are updated; the module mains still write their replays to `results/`. The note gained no lines above this one, so line references into it still hold.*

*Paths (2026-10-17, P4-9 promotion): the strip certificate (the plan, the cache records, the `layer_t14` and `centre_t14` parts, the centre identity and the merged `tc_commutator.json`) moved unchanged to `certificates/theorem_c/strip/` (`problems/P4/certificates/README.md`); `strip/` in this note means that directory.*
