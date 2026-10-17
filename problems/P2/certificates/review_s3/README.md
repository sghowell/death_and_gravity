# Independent checks for the S3 review (P2 part (iii), P_B, l = 0–19)

Reviewer's scripts and outputs for `docs/review-2026-10-17-p2-s3-adversarial.md`. Nothing here
modifies the author's files, and no script writes outside this directory.

- Python is the repository venv, `/Users/seanhowell/dev/research/death_and_gravity/.venv/bin/python`.
- Every command below was run with `nice -n 10` from this directory
  (`problems/P2/certificates/review_s3/`). The exception is pytest, which runs from the repository
  root.
- All runs used the review worktree at `db2a37eb`, the author's head.

| script | what it checks | independent of the author's code? | outputs | CPU |
|---|---|---|---|---|
| `check_certs.py` | Static consistency of all 20 `winding_PB_l*.json`. Exact tiling of the right half H: from (0, −b), via (a, −b) and (a, e), to (0, e). a = `W_l_exact` of `certificates/lg/lg_PB_l{l}.json`. b ≥ 0.55κ₋, 0 < e ≤ κ₋/10, b < κ₊. Each segment within ρ/2 of its centre. √2ρ < the distance to E₊ ∪ E_c. Exact rational re-sum of the increment balls. Also reads the Krawczyk and controls files | yes: parameters from the exact definitions in mpmath, geometry in `fractions`, no `p2` import | `check_certs.json` | 5 s |
| `indep_w.py` | Module: a float W_l from the unfactored radial equation, as the system u = R, v = ΔR′. Frobenius start for R, scipy DOP853, frozen Jost factor | yes: no `p2` import, and a formulation different from both arbw and zsol | — | — |
| `indep_winding.py 0-19` | Every stored segment re-traced with the float W and compared with its stored ball; the grid is adaptive, with steps < 0.2 rad and n/2n agreement to 1e-9. The LEFT half of ∂R′_l is traced explicitly, giving a full-boundary winding without the symmetry | yes | `indep_winding.json`, `indep_winding_stdout.txt` | 86 s |
| `replay_author.py 19,2,1,0` | The author's `part3.cauchy_segment`, unchanged, re-run with K re-derived on every stored segment of l = 0, 1, 2, 19 (306 segments). Checks overlap with the stored balls, the K and piece counts, and the re-sum | no: this checks reproducibility of the author's pipeline | `replay_author.json`, `logs/replay_l*.jsonl`, `replay_author_stdout.txt` | 2015 s |
| `recount.py 0,1,2` and `recount.py 5` | A rigorous recount around the FULL boundary, with no symmetry. The plan, the Cauchy/DFT step, the aliasing and truncation bounds, the pieces and the Arg increments are the reviewer's own. Point values come from S1's `arbw` (128 bits for l ≤ 2, 256 bits at l = 5) | partly: only zsol's box bound M is reused | `recount.json`, `logs/recount_l*.jsonl`, `recount_stdout.txt`, `recount_l5_stdout.txt` | 481 s + 1247 s |
| `box_checks.py A` | For all 2351 production discs: the float \|W\|, sampled on a 9×9 grid of the box and at 32 points of the circle, never exceeds the stored M | yes | `box_checks_A.json`, `box_checks_A_stdout.txt` | 38 s |
| `box_checks.py B` | zsol's box enclosure on 48 boxes: half-widths κ₋/16, 1e-4κ₋ and 1e-8κ₋, at high-frequency corners and interior points, for l = 1, 5, 10, 19. It must contain the float W at 16 points and S1's arbw balls at two corners | the test values are independent | `box_checks_B.json`, `box_checks_B_stdout.txt` | 73 s |
| `c0_and_krawczyk.py` | C0a: the original (`89a2336b`) C0, replayed. C0b: the test point dropped in `3e7fca05`. K1: `krawczyk_z` re-run. K2: certified sign change of W₀(iyκ₋) across the stored enclosure, with zsol and with arbw. K3: the float zero | K2 and K3 are independent of the Krawczyk code | `c0_and_krawczyk.json`, `c0_and_krawczyk_stdout.txt` | 6 s |
| `c0_detail.py` | The dropped test point under the original C0. arbw at 256 bits, with the C0 wrapper, at the three high-frequency C0 points | — | `c0_detail.json`, `c0_detail_stdout.txt` | 2 s |
| `sympy_identities.py` | The reviewer's own derivation: partial fractions of r²/Δ; the factored operator; the factored W equal to the frozen W exactly; A₊ = 1/(2κ₊) and A_c = −1/(2κ_c); the indicial roots; D4 from the Z-form; conjugation | yes | `sympy_identities.json`, `sympy_identities_stdout.txt` | < 1 s |
| `sympy_negative_controls.py` | A wrong partial-fraction coefficient, and the wrong sign of 2iωr², must both give False | yes | `sympy_negative_controls_stdout.txt` | < 1 s |
| pytest | Fast suite and P2_SLOW=1 suite | — | `replay_pytest_fast.txt`, `replay_pytest_slow.txt` | 35 s + 425 s |

Total reviewer compute: about 4450 CPU-s (1.24 core-h). At most 3 concurrent workers were used.

Replay:

```
cd problems/P2/certificates/review_s3
PY=/Users/seanhowell/dev/research/death_and_gravity/.venv/bin/python
nice -n 10 $PY check_certs.py
nice -n 10 $PY indep_winding.py 0-19
nice -n 10 $PY replay_author.py 19,2,1,0 --workers 3     # resumable (logs/replay_l*.jsonl)
nice -n 10 $PY recount.py 0,1,2 --workers 3              # resumable (logs/recount_l*.jsonl)
nice -n 10 $PY recount.py 5 --workers 3
nice -n 10 $PY box_checks.py A
nice -n 10 $PY box_checks.py B
nice -n 10 $PY c0_and_krawczyk.py
nice -n 10 $PY c0_detail.py
nice -n 10 $PY sympy_identities.py
nice -n 10 $PY sympy_negative_controls.py
cd ../../../..   # repository root
nice -n 10 $PY -m pytest problems/P2/tests -q -p no:cacheprovider                # replay_pytest_fast.txt
P2_SLOW=1 nice -n 10 $PY -m pytest problems/P2/tests -q -p no:cacheprovider     # replay_pytest_slow.txt
```

The resumable scripts skip segments already recorded in their JSONL logs. To recompute from
scratch, delete `logs/replay_l*.jsonl` or `logs/recount_l*.jsonl`.

`compare_author_jsonl.py` compares the author's gitignored working files
(`results/part3/*_segments.jsonl` in the author's worktree, read only) with the certificate pieces.
They contain 2351 ok records, one per segment, and each equals its piece
(`compare_author_jsonl_stdout.txt`). The script needs the author's worktree, since those files are
not in the repository.
