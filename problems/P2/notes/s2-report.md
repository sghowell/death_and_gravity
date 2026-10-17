# P2 S2 report: part-(i) feasibility probe at P_B

Status: S2 exit, 2026-10-17. Branch `p2/s2-part-i`, from main at `9f1b8e6f`.
**Verdict: GO for S3** under the pre-registered rule (`s2-preregistration.md` §6).
The kill criterion is not met: l = 100 certifies, and the projected l₀ is 2.

## 0. Decision record (step 0)

- Problems document revision 1.2 is approved as proposed in `s0-gate.md` §7, with
  one change: the two-sided β enclosure is a stretch goal, not part of the entry
  point.
- The .tex carries the revision-log entry and the P2 Status, Entry point and
  bibliography changes. The PDF was rebuilt with tectonic. The one new TeX
  warning is a mild underfull line in a bibliography entry.
- FORMULATION v0.3 records the decision and s0-gate.md has a one-line decision
  record (commit `3aa4af27`).
- No CLAIMS.md row depends on the decision: the rows P2-par, P2-0, P2-D1′,
  P2-2 and P2-K6 are stated at both points.

## 1. Route and design

The route is **M3**: complex rotation x = x_P + e^{iπ/4}s of the tortoise
coordinate through the photon sphere. On the rotated line a QNM becomes a decaying
solution of −u″ + e^{iπ/2}(V_l − ω²)u = 0. Exclusion follows from

  Re⟨e^{−iφ}(−u″ + Qu), u⟩ ≥ ∫ c|u′|² + q|u|² ≥ μ‖u‖² > 0,

where:
- φ is a variable phase, c = cos φ − η|φ′| and q = Re(e^{−iφ}Q) − |φ′|/(4η);
- μ comes from a ground-state (Allegretto–Piepenbrink) pointwise inequality
  E = q + c′G′ − cG′² + cG″ ≥ μ, checked in Arb on s-cells.

Proof: `notes/s2-m3.md`. At the barrier top the rotation turns the inverted
oscillator into a harmonic one, so the trapping gap γ₀ is built in: the ideal
margin is 1 − β/γ₀ = 0.673 of μ_HO = 2(l+½)Ω_Pγ₀.

The brief named M2 as the leading route and M3 as the alternate. M3 was run
first, and the choice was recorded in the pre-registration before production
(§2 there). M2's O(h)-quality centre needs the complex-coalescence Weber map of
route M1, which is open analysis; M3 needs only elementary functions. M2 was not
implemented.

## 2. Certified windows (P_B; A0–A4 all pass)

The window is 𝒲_l = [ω_l − 8γ₀, ω_l + 6γ₀] × [−β, 0], with β = 0.55κ₋. In
scaled units it runs from 8 oscillator quanta below the barrier top to 6 above.
Certificates are in `certificates/m3/`, at code commit `c4159e7f`.

| l | Re ω/κ₋ | μ_cert | μ_HO | μ_cert/μ_HO | float λ at low corner/μ_HO | min c | a_l/κ₋ ≤ top |
|---|---|---|---|---|---|---|---|
| 20 | [87.370, 110.935] | 0.32455 | 0.61436 | **0.528** | 0.661 | 0.119 | 106.64 ≤ 110.94 |
| 50 | [234.935, 258.500] | 0.85991 | 1.51342 | **0.568** | 0.681 | 0.120 | 255.51 ≤ 258.50 |
| 100 | [480.877, 504.442] | 1.68539 | 3.01185 | **0.560** | 0.679 | 0.121 | 503.43 ≤ 504.44 |

What each acceptance criterion found:

- **A0 (geometry, independent of l).**
  - The validated contour satisfies Re r*(r(s)) = x_P + s/√2 within its balls.
  - The sectors between the real axis and the contour are singularity-free. The
    proof uses the Re r* identity and Σ A_i = 0: the lower-left sector is free
    because Re r*(0) = −2.29664 > x_P = −2.43797, and the two special vertical
    segments of the upper-right sector were validated.
  - Tail bounds: |r − r₊| ≤ 0.00956 for s ≤ −64 and |r − r_c| ≤ 0.351 for s ≥ 32.
- **A1.** E ≥ μ_cert > 0 on all cells and on both tails. c ≥ 0. The Jost decay
  condition holds.
- **A2 (negative control).** Boxes of half-width γ₀/50 around each l's PS
  fundamental mode (−1.6836i, −1.6833i and −1.6832i in κ₋ units) do not
  certify: μ = −0.068, −0.097 and −0.285 with the window ground state, and the
  box-specific design fails. The method therefore cannot certify a true QNM.
  Its float margin is essentially exact: about −0.001 μ_HO at the PS mode for
  l = 100.
- **A3.** μ_cert ≤ the float optimum, for every l.
- **A4.** `python -m p2.m3 verify` replays every verdict and μ_cert identically
  (61 s).
- **P2-2 overlap.** a_l ≤ window top at every l. For l = 100 the threshold was
  computed in S2 with the S1 code (`certificates/m3/lg_PB_l100.json`, criterion
  1.881 < 2). The margin there is 1.0κ₋, which is thin but certified.

So at l = 20, 50 and 100 the whole set {Re ω ≥ ω_l − 8γ₀, −(½+δ)κ₋ ≤ Im ω ≤ 0}
is QNM-free at P_B: the window covers the part below a_l, and P2-2 the rest.

## 3. Error constant and projected l₀ (pre-registered fit)

| fit of m(l) = μ_cert/μ_HO | m_∞ | C | max residual | projected l₀ |
|---|---|---|---|---|
| m_∞ − C/(l+½) | 0.577 | 0.94 | 0.010 | **2** |
| m_∞ − C₂/(l+½)² | 0.567 | 16.0 | 0.007 | 5 |

- The certified m(l) is not monotone (0.528, 0.568, 0.560). Against the float
  eigenvalue at the window's low corner, the certificate loses about
  0.11–0.13 μ_HO to its discretization (spline ground state, mean-value cells).
- The method itself is nearly exact. At the l-th PS mode, a true QNM where the
  ideal margin is 0, the float margin is −0.020, −0.004 and −0.001 μ_HO for
  l = 20, 50 and 100. That is O(l⁻²).
- In s0-gate §4.1 terms, the barrier-top regime poses no obstruction at any
  l ≥ 20, and on this fit at any l ≥ 2. The planning range (l₀ ≈ 6–60 at P_B for
  O(1/l) methods) is beaten.

## 4. What S2 does not establish (carried to S3/S4)

- **Uniformity in l.** Only l = 20, 50 and 100 are certified. Part (i) needs the
  inequality for all l ≥ l₀ via the scaled variable σ = s(l+½)^{1/2} and closing
  boxes at h = 0. There the condition reduces to the harmonic-oscillator
  ground-state bound with margin 0.673 minus design losses; its h-structure
  matches the P4 "closing box" pattern.
- **Sub-barrier regime** (Re ω < ω_l − 8γ₀, down to the low-frequency region)
  and the **low-frequency regime** near the horizon throat. The probe did not
  address them. The float design study suggests the same contour keeps a
  positive margin at least to 8 quanta below the top. Lower energies may need a
  different phase profile, or the forbidden-region Riccati argument of s0-gate
  §4.
- **P_A** has not been probed. Its strip reaches 0.667γ₀, so the ideal M3 margin
  shrinks to 1 − 0.667 = 0.333 μ_HO. The design losses seen here (~0.12 μ_HO)
  would leave a margin of about 0.2. That is still positive, but tight.

## 5. Cost

All runs used at most 2 cores at nice 10; one P4 probe shared the machine.

| item | CPU |
|---|---|
| geometry certificate | 0.6 s |
| window certificates | 64 s (l=20, deep bisection), 0.9 s (l=50), 0.8 s (l=100) |
| negative controls | 58 / 44 / 43 s |
| a_100 | 0.9 s |
| replay | 61 s |
| tests | 12 s fast; 197 s with `P2_SLOW=1` (52 pass) |
| float design study (scratch) | ≈ 3 min |
| **S2 production total** | ≈ 4 min |

## 6. Suggested CLAIMS text (for the coordinator's bundle; not applied)

New row:

| id | problem | statement | level | evidence | date |
|---|---|---|---|---|---|
| P2-1w | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20) and for l = 20, 50, 100: W_l(ω) ≠ 0 for all ω with (l+½)Ω_P − 8γ₀ ≤ Re ω ≤ (l+½)Ω_P + 6γ₀ and −(½+δ)κ₋ ≤ Im ω ≤ 0 (the barrier-top window). Certified coercivity margins μ/μ_HO are 0.528, 0.568 and 0.560, with μ_HO = 2(l+½)Ω_Pγ₀. With P2-2 (a_l ≤ the window top; a_100 = 503.43κ₋ computed in S2) and P2-0, the region Re ω ≥ (l+½)Ω_P − 8γ₀ of S_δ is QNM-free at these l | CERTIFIED (Arb; route M3: π/4 rotation of r* through the photon sphere, variable-phase coercivity and ground-state inequality on s-cells, singularity-free sectors and tail bounds); the written implication (notes/s2-m3.md) is at WRITTEN_PROOF level, pending independent review | problems/P2/src/p2/m3.py, m3_contour.py; certificates/m3/; notes/s2-m3.md, s2-preregistration.md, s2-report.md; tests/test_m3.py | 2026-10-17 |

Optional amendment to P2-2: "also l = 100 at P_B (a_100 = 503.43κ₋,
certificates/m3/lg_PB_l100.json)".

## 7. Commits

| commit | subject |
|---|---|
| `3aa4af27` | P2: record the approval of problems-document revision 1.2 |
| `0d115788` | P2: pre-register the S2 part-(i) feasibility probe |
| `c4159e7f` | P2: implement the S2 part-(i) route M3 at the barrier top |
| (this commit) | P2: S2 barrier-top certificates, proof note and report; KAT-6 replay hygiene |

## 8. Reproduction (cwd `problems/P2/src`)

```
python -m p2.m3 run 20 50 100   # geometry, windows, controls, overlap, fit (~3 min)
python -m p2.m3 verify          # A4 replay (~1 min)
python -m pytest problems/P2/tests -q                 # from the repo root (~12 s)
P2_SLOW=1 python -m pytest problems/P2/tests -q       # ~3.5 min
```
