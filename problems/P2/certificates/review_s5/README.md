# P2 S5 review: checks and outputs

These files belong to the independent review of P2 stage S5 (the assembly P2-4 and the conversion
P2-5 at P_B), `docs/review-2026-10-26-p2-s5-adversarial.md`. They were produced on 2026-10-26 at the
review head `aaeb62dc` (branch `p2/s5-assembly`, from main `3dd0b2bc`). None of them changes a
tracked file.

Python: the repository's `.venv` (python-flint 0.9.0, mpmath 1.3.0). All runs were at `nice -n 10`.
Paths below are relative to this folder unless they say otherwise.

| file | what it is | replay (from this folder) | result |
|---|---|---|---|
| `indep_assembly.py`, `indep_assembly_out.json` | Independent re-check of the P2-4 bookkeeping. Shares no code with `p2`: exact rationals and `mpmath.iv`; the point P_B is recomputed from its definitions by exact bisection and only compared with `params.json`. Seams (b − d_δ, e, poles, γ₀, κ's), all 20 `lg`/`winding` pairs (a_l identity, metadata, criterion, path, discs, pole clearance, exact re-sum with an interval for π), l₀ at h₀ = 2/41, the regime union, a 20 295-point sampled cover of S_δ × {l} (seam points and random l up to 10⁶), β and D4 | `python indep_assembly.py` | 280/280 pass; b − d_δ = 4.54603e−10, width of the d_δ enclosure 2.7e−58 |
| `mutations.py`, `mutations_out.json` | 22 extra negative controls for the author's checker `p2.assembly.check` (beyond N1–N8), each on an in-memory copy; plus 3 rows (B1–B3) showing what the checker does not see by design | `python mutations.py` | M1–M22 all fail on the named items; B2, B3 pass (expected: analytic content is the ingredients' own replays' job); B1 is caught only incidentally |
| `quotes_check.py`, `quotes_check_out.json` | The 38 passages of [HV17] and [CF17] that the review relies on, found verbatim in the arXiv TeX, with line numbers | `python quotes_check.py --src DIR` (DIR holding `cauchyhor-r.tex`, `RN_deSitterRevised6.tex`; without `--src` it downloads the arXiv e-prints) | 38/38 found; both source sha256 match the ones recorded in the script |
| `part3_fresh_sample.py`, `part3_fresh_sample_out.json`, `part3_fresh_sample.txt` | 64 S3 winding segments (4 per l for l = 3–18) chosen off the standard deterministic sample, recomputed with the author's `cauchy_segment` and their stored K | `python part3_fresh_sample.py 2` (2 workers) | 64/64 increment balls overlap the stored ones; 350 s CPU |
| `replay_part3_PB_sample.json`, `replay_part3_PB_sample.txt` | The author's replay `p2.part3 verify` (exact re-sum and coverage of all 20 certificates, the standard 5% sample of 132 segments recomputed) | from `problems/P2`: `PYTHONPATH=src python -m p2.part3 verify --point PB --l 0-19 --sample 0.05 --workers 2 --out certificates/review_s5/replay_part3_PB_sample.json` | OK; 780 s CPU |
| `replay_params.txt`, `replay_lg_PB_0_19.txt` | `p2.params verify`; `p2.lg_bounds verify PB 19` | from `problems/P2`: `PYTHONPATH=src python -m p2.params verify`; `PYTHONPATH=src python -m p2.lg_bounds verify PB 19` | OK; OK |
| `pytest_fast.txt`, `pytest_slow.txt` | The P2 suite, fast and with `P2_SLOW=1` | from the repository root: `python -m pytest problems/P2/tests -q -p no:cacheprovider` (and with `P2_SLOW=1`) | 90 passed, 17 skipped (37.5 s); 107 passed, 1 pre-existing warning (685 s) |

**Follow-up, 2026-10-26** (the author's fix commit `3f2a17aa`; review §9):

| file | what it is | replay | result |
|---|---|---|---|
| `followup_check.py`, `followup_check_out.json`, `followup_tex_item_v1.3.txt` | Mechanical checks of the fix commit: the s5-report §6 rows are byte-identical to the review's §7 rows; the .tex `\item[v1.3]` is the pre-review text plus exactly the review's (c), (d) corrections and the two dates; the pre-registered plan is a byte prefix of the current plan; no residual scope wording; no code, test, certificate, `docs/` or `CLAIMS.md` change. Writes the final `\item[v1.3]` with the one-word change of §9 (5(i)) | `python followup_check.py` | all pass |
| `followup_pytest_fast.txt`, `followup_assembly_verify.txt` | The fast P2 suite and `p2.assembly verify` at `3f2a17aa` | as above (fast suite; `assembly verify`) | 90 passed, 17 skipped; identical |

The assembly replay itself, `PYTHONPATH=src python -m p2.assembly verify` from `problems/P2`, prints
`verify: identical` (2.8 s); it writes nothing.

Sources read by the review (fetched from arXiv as TeX on 2026-10-26):
- [HV17] arXiv:1512.08004v2 (posted 2020-05-27; TeX date line "December 25, 2015. Final revision: July
  12, 2017."), `cauchyhor-r.tex`, sha256 `a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be`.
  J. Math. Phys. 58, 081509 (2017).
- [CF17] arXiv:1607.01018v2 (posted 2017-10-23), `RN_deSitterRevised6.tex`, sha256
  `dcaab0e03e43622592f6f2d336bf131a36085d627d62efdf252903a37464d453`. Ann. Henri Poincaré 18,
  3371–3398 (2017), DOI 10.1007/s00023-017-0592-z.
