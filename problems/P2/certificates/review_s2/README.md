# Reviewer scripts for the P2 S2 review (2026-10-17)

These files back `docs/review-2026-10-17-p2-s2-adversarial.md`. They were written by the reviewer, not the
S2 author. They read the author's certificates and, where a comparison needs it, call the author's code.
The independent checks do not reuse the author's code paths: contour, potential, derivatives, cells and
tails are all the reviewer's own. Python is the repository venv (`.venv/bin/python`), run from the repo
root. Captured stdout is in `*_stdout.txt` and `replay_*.txt`; machine-readable results are in the
matching `*.json`.

| script | what it checks | method (independent of the author?) | run time |
|---|---|---|---|
| `indep_geometry.py` | **Item 1a/1b.** Parameters; x_P < Re r\*(0); the contour image r(s), s ∈ [−64, 32]; the two special vertical segments and a fan of vertical lines over both sectors; the nearest r = 0 and r = ∞ points; the tail discs | mpmath parameters checked against `params.json`. Krawczyk enclosures of the branch-tracked r\*(r) = x_P + e^{iπ/4}s, compared with the author's Taylor-model balls at all 3,074 points of spacing 1/32. scipy DOP853 for the vertical lines. Triangle-inequality tail discs in Arb (yes, except that the author's `Ray` is called to get the balls being compared) | 1.4 s |
| `indep_identities_sympy.py` | **Items 1c/1d and code formulas.** The integration-by-parts identity of §3, Young's step, the Lemma 4 ground-state identity, and E′, c′, c″, q′, Q′, V, V′ and F as coded in `m3.py`. Also the rotated equation, the decay exponent (a+b)/√2, the hull monotonicity, the tail monotonicity of \|φ′\| and the barrier-top model | sympy (yes) | 0.4 s |
| `indep_cells_arb.py` | **Item 2.** E ≥ μ_cert on all of s ∈ [−64, 32] and both tails, for every ω² in the reviewer's own hull of the window, with the stored ground-state nodes; c > 0; decay condition | Arb 128-bit: `params.json` balls, Krawczyk r(s) (point and interval), expanded-polynomial V_l, dual-number derivatives, own cover and bisection, subdivided tail discs (yes) | 1.5 s |
| `indep_float_spectrum.py` | **Items 1a, 3, 5.** (A) complex scaling: the eigenvalues of −d²/ds² + iV_l on the rotated line reproduce the PS modes, and none lies in a window. (B) the float ground-state margin at the four window corners, along the bottom edge, on the column Re ω = ω_l and at the PS mode. (C) the float min of E with the stored ground state over the ω² hull and over the true window corners | numpy/scipy FD (yes) | 3.8 s |
| `sharp_negative_control.py` | **Item 3.** The author's rigorous checker with a ground state designed near the PS mode (Im ω = −0.88γ₀), on boxes that contain / straddle / avoid the mode | the author's `m3.certify_l` (deliberately: it tests that code path) | 416 s for l = 100, 50, 20 (deep bisection) |
| `indep_p22_overlap.py` | **Item 4.** a_l values equal the lg certificates; certified a_l < window top with the reviewer's window; the S1 criterion replayed at l = 100; the actual \|W e^{iΘ}/(2iω) − 1\| at l = 100 from the reviewer's own real-axis ODE integration | Arb comparisons on `params.json` balls; `lg_bounds.LG` for the replay; scipy DOP853 (yes, except the replay) | 28 s |
| `float_leak_tail_bound.py` | **Item 2.** The size of the float round trip in `m3.tail_bound` (m3.py:398–399), and the tail checks redone with the disc enlarged | the author's routine against an exact-Arb version | 1 s |

Replays of the author's commands:
- `replay_m3_verify.txt`: `python -m p2.m3 verify`, cwd `problems/P2/src`.
- `replay_pytest_fast.txt`: `pytest problems/P2/tests`.
- `replay_pytest_slow.txt`: the same with `P2_SLOW=1`.

All were run at nice 10, and `git status` showed no change to tracked files afterwards.

```
.venv/bin/python problems/P2/certificates/review_s2/indep_geometry.py
.venv/bin/python problems/P2/certificates/review_s2/indep_identities_sympy.py
.venv/bin/python problems/P2/certificates/review_s2/indep_cells_arb.py          # [l ...], default 20 50 100
.venv/bin/python problems/P2/certificates/review_s2/indep_float_spectrum.py
.venv/bin/python problems/P2/certificates/review_s2/indep_p22_overlap.py
.venv/bin/python problems/P2/certificates/review_s2/float_leak_tail_bound.py
.venv/bin/python problems/P2/certificates/review_s2/sharp_negative_control.py   # [l ...]; the run used 100 50 20
```
