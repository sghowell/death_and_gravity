# P2 S6 review: checks and outputs

These files belong to the independent review of P2 stage S6 (the two-sided β at P_B, P2-6, and
Lemma E, P2-E), `docs/review-2026-10-27-p2-s6-adversarial.md`. They were produced on 2026-10-27 at the
review head `172eb48c` (branch `p2/s6-beta`, from main `300c22a7`), in the worktree
`.claude/worktrees/p2-s6-review` (branch `p2/s6-review`). None of them changes a tracked file.

Python: the repository's `.venv` (python-flint 0.9.0, mpmath 1.3.0, sympy 1.14.0, numpy 2.4.6).
All runs were at `nice -n 10`, with at most three of the reviewer's processes at once. Paths are
relative to this folder unless they say otherwise. The scripts of the first table import nothing from
`p2` (they read certificates only as data), except `mutations.py`, which drives the author's checker.

## Independent checks (reviewer code)

| file | what it is | replay (from this folder) | result |
|---|---|---|---|
| `indep_rnds.py` | The reviewer's toolkit (python-flint Arb only): P_B from its definitions (interval Newton for r_ext and the four roots of Δ), the phase-factored radial equations, validated Frobenius series at r₊ and r_c and validated Taylor continuation, each with a geometric-majorant tail, and W_l in the form of FORMULATION §1.2 | (library) | — |
| `indep_lemmaE.py`, `indep_lemmaE_out.json`, `indep_lemmaE_log.txt` | Lemma E, rigorously (200 bits): W̃_l(−iκ₊) for l = 0, 1, 2, 3, 5, 10, 19, 20, 50, 100 (exponent-1 solution Z2 with B₀ = 0 exactly); the residue check ε·W_l(ω₀ + ε) against i(2κ₊r₊ + L)/(2r₊²)·W̃_l at ε = 10⁻³, 10⁻⁵, 10⁻⁷ κ₋e^{iπ/4} (l = 0, 1, 5, 19); a sign change of the real function W₀(−iyκ₋) across P2-K6's ball; P_B against `params.json` | `python indep_lemmaE.py` | pass. W̃₀ = 0.0382855818144094394 ± 2·10⁻²²; every W̃_l overlaps the author's ball; signs + (l = 0), − (l ≥ 1); the residue error is O(ε) (ratios 99.7–100.0 per factor 100 in ε); W₀ changes sign (−1.50·10⁻³⁰ → +4.49·10⁻³⁰) across the ball |
| `sympy_lemmaE.py`, `sympy_lemmaE_out.json` | Lemma E's algebra in sympy: the phase-factored equations from (R); D_m = Δ₁m(m − iω/κ₊); the first step; the log coefficient (2κ₊r₊ + L)/Δ₁ of the exponent-0 solution at ω₀; the residue limit; the Wronskian form; the identities behind the closed form (E-int); the limits of the Wronskians at r_c; W6′'s scaled log coefficient ã₁(ω₀) = h²(2κ₊r₊ + L)/Δ′(r₊) (E9) | `python sympy_lemmaE.py` | all true |
| `indep_float.py`, `indep_float_out.json`, `indep_float_log.txt` | Float argument-principle counts with the reviewer's own W_l (numpy, vectorised; adaptive steps), checked against `indep_rnds` at 1500 bits to ≤ 6·10⁻¹⁴ relative at 9 points up to Re ω = 131κ₋: the full boundary and the right half of every band box B_l (l ≤ 19); small boxes around −iκ₊ (l = 0, 1, 7, 19) and ω_*; ω_* by Newton; the l = 20 and l = 21 boxes [−1.3(l+½)Ω_P, 1.3(l+½)Ω_P] × [−b′, −b₀] (reaching part (i)'s regime A) | `python indep_float.py` | pass (51 s). Windings 0 (l = 0) and −1 (l = 1–19), right half/π equal; −1 around −iκ₊; +1 around ω_*; 0 around both; −1 at l = 20, 21 (the pole, no zero); ω_* = −0.909380978994807iκ₋ |
| `indep_assembly.py`, `indep_assembly_out.json` | The seams and the cover of S_e (Arb 256 bits, exact rationals): E-points in S_e; b′ and b₀ (least/greatest 2⁻³⁰ dyadics); a′_l against the A certificates and ζ_AΩ_P/h_l; ω₀ interior; the traced paths, discs and pole clearance against every n of E₊ ∪ E_c; exact re-sums; ω_* in B₀; W7 at each fixed h_l; part (i)'s cover inequalities and tilings; a 40 000-point sampled cover of S_e for l up to 10⁶ | `python indep_assembly.py` | 215/215 pass; b′ − eκ₋ = 3.879·10⁻¹⁰ (eκ₋ known to ~10⁻⁶⁰); min pole clearance 1.517·√2ρ; ω_* 0.010619κ₋ above B₀'s bottom; no sample uncovered |
| `mutations.py`, `mutations_out.json` | 36 further negative controls for `p2.beta_assembly.check` (in-memory copies), and 3 rows (B1–B3) of what the check does not see by design | `python mutations.py` | M1–M36 all caught on the named items; B1 (an increment within the ¼ tolerance) and B2 (a `+dirty` tag) pass, as designed; B3 (a disc far from E inflated by 5/4) is rejected only by the conservative n ≥ 3 bound 3 min(κ₊, κ_c) + Im c > √2ρ of S5's `_poles_ok`, which ignores Re c (the stored pieces pass that bound with margin ≥ 1.095) |
| `quotes_check.py`, `quotes_check_out.json` | 16 passages of [HV17] v2 and [Vas13] v2 used for the notion of resonance at −iκ₊, found verbatim in the arXiv TeX with line numbers | `python quotes_check.py --src DIR` (DIR holding `cauchyhor-r.tex` and `kerr-ah-b.tex`; without `--src` it downloads the e-prints) | 16/16 found; both source sha256 match |
| `dirty_episode_check.py`, `dirty_episode_check_out.json` | The `+dirty` episode: the author's 178 first-pass L rows (read, not written, from the git-ignored `results/beta/part1/` of the author's worktree) against the committed `L_PB.json` | `python dirty_episode_check.py` | 178/178 leaves identical; the final 1280 rows all tagged `7168b2e4` clean |

## Replays of the author's certificates (author's code, outputs written only here)

`replay_author.py` calls the author's library functions directly (not the run/launch drivers), so
nothing under `certificates/beta` or `results/` is written. Run it from this folder; it puts
`problems/P2/src` on `sys.path` itself.

| job | what | command | result |
|---|---|---|---|
| `replay_epole.json` | Lemma E certificate recomputed | `python replay_author.py epole` | identical apart from timing and commit |
| `replay_part2.json` | part (ii′), A at fixed h_l, all 20 l | `python replay_author.py part2` | identical (66 s) |
| `replay_TA.json` | part (i) T (10 boxes), A, and the cover | `python replay_author.py TA` | identical |
| `replay_controls.json` | NB3–NB5 | `python replay_author.py controls` | identical |
| `replay_L.json` | all 1280 L grid boxes (pole-free series) | `python replay_author.py L` | 1280/1280 leaves identical to `L_PB.json` (243 s) |
| `replay_S_plus.json`, `replay_S_c.json` | all S jobs, both sides (extended grid) | `python replay_author.py S+`, `... Sc` | 45/45 leaves exactly identical on each side, float margins included (197 s, 34 s) |
| `replay_band.json` | 72 band pieces: the 8 nearest to −iκ₊ and ω_* (l = 0, 1) and 64 drawn with seed 20261027 from the pieces outside the author's 5% replay; recomputed with `part3.cauchy_segment` at the stored K | `python replay_author.py band 20261027 64` (2 workers) | 72/72 overlap (377 s CPU) |
| `replay_beta_assembly_verify.txt`, `replay_p2_assembly_verify.txt` | the author's verifies | from `problems/P2`: `PYTHONPATH=src python -m p2.beta_assembly verify`; `... -m p2.assembly verify` | identical; identical |
| `pytest_fast.txt`, `pytest_slow.txt` | the P2 suite, fast and with `P2_SLOW=1` | from the repository root: `python -m pytest problems/P2/tests -q -p no:cacheprovider` (and with `P2_SLOW=1`) | fast: 102 passed, 21 skipped (67 s); slow: 123 passed, 1 pre-existing warning (S2's `m3` overflow), 732 s |

Each `*_log.txt` holds the run's stdout and `/usr/bin/time -p`.

## Sources read by the review (fetched from arXiv as TeX on 2026-10-27)

- [HV17] arXiv:1512.08004v2 (posted 2020-05-27; TeX date line "December 25, 2015. Final revision: July 12,
  2017."), `cauchyhor-r.tex`, sha256 `a0d3302f755222eb493170e1eea60f988a44c96f0e69108c89171e99733cf2be`.
- [Vas13] arXiv:1012.4391v2 (updated 2011-11-09; TeX date line "Revised: November 1, 2011. Original version:
  December 20, 2010"), `kerr-ah-b.tex`, sha256 `b7ead254e1b5aa3654d8ba184a89f5fd4a1fb6a7e2cc79c00356f3030b83bc18`.
