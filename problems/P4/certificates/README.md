# P4 certificates

The machine-checkable artifacts behind P4's CERTIFIED ledger rows, as the main
`README.md` requires ("Ledger levels": a CERTIFIED row's certificate is
committed under `certificates/`).

Until 2026-10-17 P4 kept these files under `problems/P4/results/`. They were
moved here with `git mv`, unchanged. `SHA256SUMS` lists every file; check them
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
| `theorem_c/review_s45/` | P4-8 (and P4-9) | the S4-5 review's independent checks (its `README.md`) | the review, `docs/review-2026-10-17-p4-s4-5-adversarial.md` |
| `recert_2026-10-17/` | P4-3, P4-6, P4-8 | re-runs of all of the above on the fixed code (`eba80cb4`, `13119db1`), and `comparison.json` | `scripts/p4_recertify.py run` / `compare` |

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

These sources are hash-pinned by the 1092 `layer_t14` parts behind P4-9, so
they are not edited. A replay therefore writes a fresh file under `results/`
and never overwrites a certificate.

`scripts/p4_recertify.py run OUTDIR` re-runs everything into `OUTDIR`, and
`compare OUTDIR` checks the result against the files here.

`scripts/s43_parts.py` writes `checks.json`, `lemma_d_certified.json` and
`p43_centre.json` here directly, and its merge reads `checks.json` from here.

## Not here yet

P4-9 (strip exclusion) and P4-10 (Theorem C) are CONJECTURED until the 30
centre pieces land. Their artifacts stay under `results/theorem_c/` until the
rows are promoted:

- `parts/layer_t14/`;
- `layer_recompute/merge_substitute_layer_t14.json`;
- the S4-3 plan and caches.

## Old paths (2026-10-17 move)

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

Dated records keep the old paths: the review documents in `docs/` and the
paths stored inside these JSON files. The P4 notes and the ledger use the new
ones.
