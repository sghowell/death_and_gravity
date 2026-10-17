# P4 certificates

The machine-checkable artifacts behind P4's CERTIFIED ledger rows, as the main
`README.md` requires ("Ledger levels": a CERTIFIED row's certificate is
committed under `certificates/`).

Until 2026-10-17 P4 kept these files under `problems/P4/results/`. They were
moved here with `git mv`, unchanged: first the files of P4-3, P4-6 and P4-8,
then the strip certificate when P4-9 and P4-10 were promoted (`theorem_c/strip/`,
`theorem_c/review_s46/`). `SHA256SUMS` lists every file; check them
with `cd problems/P4/certificates && shasum -a 256 -c SHA256SUMS`.

## Files and ledger rows

| file | row | what it certifies | produced by |
|---|---|---|---|
| `theorem_b/winding_R.json` | P4-6 | winding number 2 of E on ∂R, 290 contour tiles | `python -m p4.validated.modecount TUBE_B` |
| `theorem_b/analyticity_R.json` | P4-6 | analyticity cover of R, 432 boxes | `python -m p4.validated.analyticity TUBE_B` |
| `theorem_b/krawczyk.json` | P4-6 | Krawczyk enclosures of κ̄ and κ₁, and γ = 1/κ₁ | snippet in `notes/s2-theorem-b.md` §4 |
| `theorem_c/tc_monotone.json` | P4-3 | tube, sonic and centre-tail sign conditions | `python -m p4.validated.tc_monotone TUBE` |
| `theorem_c/p43_centre.json` | P4-3 | the centre tail on the certified (a\*, μ\*) balls, with provenance | `scripts/s43_parts.py p43-centre` |
| `theorem_c/tc_energy.json` | P4-8 | tube 14520 + sonic 16666 + closing box on x ∈ [−4.5, 0] | `python -m p4.validated.tc_energy all TUBE 6.0` |
| `theorem_c/lemma_d_certified.json` | P4-8 | Lemma D on y ∈ (0, YCUT], 1100/1100 subranges, certified balls | `scripts/s43_parts.py lemma-d` |
| `theorem_c/checks.json` | P4-8, P4-9 | exact identities E1/E2/E3/E5 | `scripts/s43_parts.py checks` |
| `theorem_c/review_s45/` | P4-8, P4-9 | the S4-5 review's independent checks (its `README.md`) | the review, `docs/review-2026-10-17-p4-s4-5-adversarial.md` |
| `theorem_c/strip/plan.json` | P4-9 | every piece of the strip cover (576 pieces) | `scripts/s43_parts.py plan` |
| `theorem_c/strip/caches/{full,gap,deep,band}.json` | P4-9 | build records of the tube caches, with their content hashes | `scripts/s43_parts.py build-caches` |
| `theorem_c/strip/layer_t14/` | P4-9 | 546 S4-4 layer-engine parts on x ∈ [−4.56, 0], τ ≥ 14, ε ∈ [0, 6] | `scripts/s44_recompute.py launch --workers 4 --nice 10` (about 3 h 52 min on 4 workers) |
| `theorem_c/strip/centre_t14/` | P4-9 | 30 S4-6 centre parts on y ∈ (0, YCUT], τ ≥ 14, ε ∈ [0, 6] | `scripts/s46_centre.py launch --workers 4 --nice 10 --tau-mins 224,14 --max-evals 150000` (about 0.7 core-h) |
| `theorem_c/strip/centre_nb_identity.json` | P4-9 | the exact centre layer identity behind the raise of nB's offsets | `scripts/s46_centre.py identity --out FILE` |
| `theorem_c/strip/tc_commutator.json` | P4-9, P4-10 | the merge: status complete, `strip_exclusion_certified = true`, 576 pieces, exact cover | `scripts/s43_parts.py merge --substitute-layer-t14 --substitute-centre` |
| `theorem_c/review_s46/` | P4-9 | the S4-6 centre review's independent checks (its `README.md`) | the review, `docs/review-2026-10-17-p4-s4-6-centre-adversarial.md` |
| `recert_2026-10-17/` | P4-3, P4-6, P4-8 | re-runs of the P4-3, P4-6 and P4-8 files on the fixed code (`eba80cb4`, `13119db1`), and `comparison.json` | `scripts/p4_recertify.py run` / `compare` |

P4-10 (Theorem C) adds no file of its own: it combines P4-6, P4-8 and P4-9
through Lemma G1.

`TUBE` is the pinned cache `data/cache/tube_full.json`, and `TUBE_B` the
Theorem B tube (both gitignored; content hashes in the records). The `lemma_d`
block inside `theorem_c/tc_energy.json` is the 2026-08-31 midpoint record, which
does not enclose the EC background. Cite `theorem_c/lemma_d_certified.json`
instead. That block is also the only difference that `recert_2026-10-17/comparison.json`
reports (`tc_energy: identical = false`); every other certified value, margin
and count re-ran identically.

P4-4 and P4-5 (Theorem A) store no artifact: their certificates are the exact
and Arb tests (`tests/test_validated_*.py`, `tests/test_derivation.py`) cited
in the ledger.

## Replays

The module mains in `src/p4/validated` still write to `results/`:

- `tc_energy.py` and `tc_monotone.py` write to `results/theorem_c/`;
- `modecount.py`'s default `--out` is under `results/theorem_b/`.

These sources are hash-pinned by the 1092 `layer`/`layer_t14` parts and the 60
centre parts, so they are not edited. Nor are the pinned drivers
`scripts/s44_layer.py` and `scripts/s46_centre.py`: they still write new parts
under `results/theorem_c/parts/` and the identity under `results/theorem_c/`
(or wherever `identity --out` says). A replay therefore writes a fresh file
under `results/` and never overwrites a certificate; compare it with the file
here.

`scripts/p4_recertify.py run OUTDIR` re-runs the P4-3, P4-6 and P4-8 files into
`OUTDIR`, and `compare OUTDIR` checks the result against the files here.

`scripts/s43_parts.py` is not pinned. It writes `checks.json`,
`lemma_d_certified.json` and `p43_centre.json` here directly, and `plan`,
`build-caches` and a complete `merge` write into `theorem_c/strip/`. Its merge
reads `checks.json`, the plan, the cache records, the `layer_t14` and
`centre_t14` parts and the identity from here (`STRIP`, `layer_dir_path`).

## The strip certificate (P4-9 and P4-10, promoted 2026-10-17)

- The S4-6 centre review found no blocking or should-fix flaw and recommended
  CERTIFIED for both rows, with this move in the same commit.
- A re-run of the merge at the new paths reproduces `strip/tc_commutator.json`
  except `provenance` (git head, time, driver hash, the new `strip_dir`) and
  the identity's recorded path. The committed file is the reviewed one
  (`6c228a73…`).
- The aggregate command of `notes/s4-5-assembly.md` §8 hashes paths as well
  as contents. At the new paths it gives `f21d32bd…` for `strip/layer_t14/`
  and `7b087f9a…` for `strip/centre_t14/`. The aggregates quoted in the notes
  and reviews, `e9835a41…` and `4ac8a7b0…`, are the same files at their old
  paths: map the paths back before the outer `shasum` and they reproduce.
- Still under `results/theorem_c/`, supporting and not cited for the theorem:
  the second τ modes (`parts/layer/`, `parts/centre_t224/`), the stale S4-3
  parts, `layer_recompute/` (the 546-piece merge of 2026-10-17, historical),
  `centre_sampling.json` and `centre_probe/`.

## Old paths

| old path (`problems/P4/...`) | new path (`problems/P4/certificates/...`) |
|---|---|
| `results/theorem_b/{winding_R,analyticity_R,krawczyk}.json` | `theorem_b/` (same names) |
| `results/theorem_c/tc_energy.json` | `theorem_c/tc_energy.json` |
| `results/theorem_c/tc_monotone.json` | `theorem_c/tc_monotone.json` |
| `results/theorem_c/parts/checks.json` | `theorem_c/checks.json` |
| `results/theorem_c/parts/lemma_d_certified.json` | `theorem_c/lemma_d_certified.json` |
| `results/theorem_c/parts/p43_centre.json` | `theorem_c/p43_centre.json` |
| `results/theorem_c/review_s45/` | `theorem_c/review_s45/` |
| `results/recert_2026-10-17/` | `recert_2026-10-17/` |
| `results/theorem_c/parts/plan.json` | `theorem_c/strip/plan.json` (P4-9 promotion) |
| `results/theorem_c/parts/caches/` | `theorem_c/strip/caches/` |
| `results/theorem_c/parts/layer_t14/` | `theorem_c/strip/layer_t14/` |
| `results/theorem_c/parts/centre_t14/` | `theorem_c/strip/centre_t14/` |
| `results/theorem_c/centre_nb_identity.json` | `theorem_c/strip/centre_nb_identity.json` |
| `results/theorem_c/tc_commutator.json` | `theorem_c/strip/tc_commutator.json` |
| `results/theorem_c/review_s46/` | `theorem_c/review_s46/` |

Dated records keep the old paths: the review documents in `docs/`, the
review checks' own scripts and READMEs, and the paths stored inside these JSON
files. Replay a review's scripts on a checkout from before its move (S4-6:
`2bda78ba`), or map the paths. The P4 notes and the ledger use the new ones.
