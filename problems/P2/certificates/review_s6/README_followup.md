# P2 S6 review follow-up (2026-10-27): checks and outputs

These files belong to the dated follow-up (§9) of `docs/review-2026-10-27-p2-s6-adversarial.md`. They were
produced at the author's fix commit `43eed599` (fast-forwarded onto `p2/s6-review`), at `nice -n 10`, with the
repository's `.venv`. `README.md` in this folder describes the main review's files and is left as committed.
Paths are relative to this folder.

| file | what it is | replay (from this folder) | result |
|---|---|---|---|
| `followup_check.py`, `followup_check_out.json` | F1 the commit's scope (only FORMULATION.md and the three S6 notes; CLAIMS.md and the .tex unchanged since main); F2 the s6-report §9 rows and notes against the review's §8 (byte for byte); F3 FORMULATION: the pre-fix file with only the v1.1 pronoun change, plus v1.2 equal to §8's text with the date and the interim wording; F4 the pre-registration's old text is a byte prefix (erratum appended only); F5 the line numbers cited in s6-proof §2.2 against `quotes_check_out.json`; F6 the v1.4 draft's numbers and LaTeX; F7 writes the final v1.4 text; F8 compiles the .tex as it is, with the draft and with the final text (tectonic, in a temporary directory outside the repository) | `python followup_check.py` | all pass. The draft adds one 92.2 pt overfull line; the final text adds none |
| `followup_tex_item_v1.4.txt` | the final `\item[v1.4]` (the author's draft, with one wording fix in (c) and a typesetting-only `{\emergencystretch=3em … \par}` wrapper), leading ⟨Month day⟩ left for the coordinator | (written by `followup_check.py`) | — |
| `followup_pytest_fast.txt`, `followup_beta_assembly_verify.txt`, `followup_p2_assembly_verify.txt` | the fast P2 suite and both assembly verifies at `43eed599` | as in `README.md` (from `problems/P2` for the verifies, from the repository root for pytest) | 102 passed, 21 skipped (65 s); identical; identical |
