# Independent checks for the S3 review (P2 part (iii), P_B, l = 0–19)

Reviewer's scripts and outputs for `docs/review-2026-10-17-p2-s3-adversarial.md`. Nothing here
modifies the author's files. Python is the repository venv
(`/Users/seanhowell/dev/research/death_and_gravity/.venv/bin/python`); every command below was run
with `nice -n 10` from this directory (`problems/P2/certificates/review_s3/`), except pytest, which
runs from the repository root.

| script | what it checks | independent of the author's code? | output | CPU |
|---|---|---|---|---|
| `check_certs.py` | static consistency of all 20 `winding_PB_l*.json`: exact tiling of the right half H from (0, −b) via (a, −b), (a, e) to (0, e); a = `W_l_exact` of `certificates/lg/lg_PB_l{l}.json`; b ≥ 0.55κ₋, 0 < e ≤ κ₋/10, b < κ₊; each segment within ρ/2 of its centre; √2ρ < distance to E₊ ∪ E_c; exact rational re-sum of the increment balls; the Krawczyk and controls files | yes (mpmath parameters from the exact definitions, `fractions` geometry; no `p2` import) | `check_certs.json` | 5 s |
| `indep_w.py` | module: a float W_l from the unfactored radial equation as the system u = R, v = ΔR′ (Frobenius start for R, scipy DOP853, frozen Jost factor) | yes (no `p2` import; different formulation from both arbw and zsol) | — | — |
| `indep_winding.py 0-19` | every stored segment re-traced with the float W (adaptive grid, steps < 0.2 rad, n/2n agreement to 1e-9) and compared with its stored ball; the LEFT half of ∂R′_l traced explicitly; full-boundary winding without the symmetry | yes | `indep_winding.json`, `indep_winding_stdout.txt` | 86 s |
| `replay_author.py 19,2,1,0` | the author's `part3.cauchy_segment`, unchanged, re-run on every stored segment of l = 0, 1, 2, 19 (306 segments) with K re-derived; overlap with the stored balls; K and piece counts; re-sum | no (reproducibility of the author's pipeline) | `replay_author.json`, `logs/replay_l*.jsonl`, `replay_author_stdout.txt` | see review §6 |
| `recount.py 0,1,2` | rigorous recount around the FULL boundary (no symmetry) with the reviewer's own plan, Cauchy/DFT, aliasing and truncation bounds, pieces and Arg increments; point values from S1's `arbw` | partly: only zsol's box bound M is reused | `recount.json`, `logs/recount_l*.jsonl` | see review §6 |
| `box_checks.py A` | for all 2351 production discs: the float |W| sampled on a 9×9 grid of the box and 32 points of the circle never exceeds the stored M | yes | `box_checks_A.json` | see review §6 |
| `box_checks.py B` | zsol's box enclosure on boxes of half-width κ₋/16, 1e-4κ₋, 1e-8κ₋ at high-frequency corners (l = 1, 5, 10, 19) contains the float W at 16 points and S1's arbw balls at two corners | the test values are independent | `box_checks_B.json` | see review §6 |
| `c0_and_krawczyk.py` | C0a: the original (89a2336b) C0 replayed; C0b: the test point dropped in 3e7fca05; K1: `krawczyk_z` re-run; K2: certified sign change of W₀(iyκ₋) across the stored enclosure (zsol and arbw); K3: float zero | K2 and K3 are independent of the Krawczyk code | `c0_and_krawczyk.json` | see review §6 |
| `sympy_identities.py` | the reviewer's own derivation: partial fractions of r²/Δ; the factored operator; the factored W equals the frozen W exactly; A₊ = 1/(2κ₊), A_c = −1/(2κ_c); indicial roots; D4 from the Z-form; conjugation | yes | `sympy_identities.json` | < 10 s |

Replay:

```
cd problems/P2/certificates/review_s3
PY=/Users/seanhowell/dev/research/death_and_gravity/.venv/bin/python
nice -n 10 $PY check_certs.py
nice -n 10 $PY indep_winding.py 0-19
nice -n 10 $PY replay_author.py 19,2,1,0 --workers 3     # resumable (logs/replay_l*.jsonl)
nice -n 10 $PY recount.py 0,1,2 --workers 3              # resumable (logs/recount_l*.jsonl)
nice -n 10 $PY box_checks.py A
nice -n 10 $PY box_checks.py B
nice -n 10 $PY c0_and_krawczyk.py
nice -n 10 $PY sympy_identities.py
cd ../../../..   # repository root
nice -n 10 $PY -m pytest problems/P2/tests -q -p no:cacheprovider
P2_SLOW=1 nice -n 10 $PY -m pytest problems/P2/tests -q -p no:cacheprovider   # replay_pytest_slow.txt
```

The resumable scripts skip segments already in their JSONL logs; delete `logs/replay_l*.jsonl` or
`logs/recount_l*.jsonl` to recompute from scratch. None of the scripts writes outside this directory.
