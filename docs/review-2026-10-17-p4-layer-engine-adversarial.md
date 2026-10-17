# Independent adversarial review: P4 S4-4 layer engine and the S4-3 substitution

Date: 2026-10-17. Reviewer: a separate agent with no part in S4-3 or S4-4.
Reviewed: branch `p4/theorem-c-s4-3` at `c3beffc2`. HEAD `c4171509` adds
parts only.

How the review ran:
- Read-only. Computations ran in a scratch copy (`git archive c3beffc2`)
  with `PYTHONPATH` set to the copy, at nice 10, one process at a time, about
  0.25 core-h in all.
- The worktree changed only through files the running workers wrote.

**Overall verdict.** The substitution is sound for the strip exclusion.
- On its piece × ε ∈ [0, 6] × τ ∈ [14, ∞), a passing `layer_t14` part
  certifies exactly what the strip argument uses: dT, dB, nB ≠ 0 and S < 0.
  That is the per-κ argument of `s4-commutator.md` §3, together with its
  endgame (the `s4-energy.md` §1 template, per-κ by its statement at
  :78–79).
- Every polynomial the endgame needs nonzero at a strip κ divides
  dT = Dp·SDr·D2r·Bcr, up to powers of Dp. This was checked exactly.
- `analyze_piece`'s per-piece, κ-uniform bounds (Rbar < 14 and
  |A|, |F_N| < 14) are not needed. In the strip the real-pole bound is
  vacuous: D2r's roots A and −F_N are real. Rbar is only a sufficient
  condition for Bcr ≠ 0, which dT ≠ 0 gives directly. Outside the strip,
  nothing uses the S4-3 pieces.
- So the literal claim of §13.2 ("certifies everything `analyze_piece`
  certifies") is false, but only as wording.

The enclosure of Σ is correct:
- the §3 identity (sympy residual 0);
- the top-slot drops and the (tt·dT_k)′ raise;
- the pivot forms;
- the cover tiling, including both pseudo-boxes and the closing piece;
- Arb's rejection of zero-containing denominators, including inside mv3;
- the mv3 mean-value argument, against the TM semantics `tc_tm` actually
  implements.

Two latent defects were found; neither was triggered (336,196 instrumented
strip evaluations, 0 misses):
- a formal flaw in shared `tc_tm` code: the evaluation ball need not contain
  the strip endpoints;
- a robustness gap in mv3: the ι-centre is computed in floats.

Independent checks:
- 17,659 independent point values on 7 pieces never exceed a certified bound;
- the bounds are 0.50–0.95 of the float suprema;
- four pieces recompute bit for bit in both modes.

With the must-fix items of §5 done, the parts merit CERTIFIED, the same level
the two-phase S4-3 path would merit.

## 1. Verdicts

| # | Claim | Verdict |
|---|---|---|
| 1 | A passing `layer_t14` part implies everything the strip argument uses on its piece | PROVEN AS STATED |
| 1a | §13.2: "certifies everything `analyze_piece` certifies … plus Rbar and real-pole checks" | GAP (wording only) |
| 1b | η-elimination: Bcr ≠ 0 at strip κ | PROVEN AS STATED |
| 1c | Real poles stay off the strip | PROVEN AS STATED (vacuous for τ ≥ 14) |
| 1d | Endgame nonvanishing (Rbar and real-pole bounds not needed downstream) | PROVEN AS STATED |
| 1e | dT, dB and nB are literally the `build_pairs` polynomials, sonic raises included | PROVEN AS STATED |
| 1f | ε and τ ranges; boundary terms; function space | PROVEN AS STATED (κ-independent, unchanged) |
| 2a | §3 normalization identity | PROVEN AS STATED |
| 2b | Top-slot drops; the (tt·dT_k)′ raise | PROVEN AS STATED |
| 2c | Pivot forms; negative ℓ-powers | PROVEN AS STATED (one sentence of the note is wrong) |
| 2d | Box geometry; the cover tiles the domain | PROVEN AS STATED |
| 2e | A zero-containing denominator rejects the box, mv3 included | PROVEN AS STATED |
| 2f | Outward rounding | FLAW (formal, latent, shared with S4-3) |
| 3a | mv3 argument and implementation | PROVEN AS STATED |
| 3b | r = Σ rad(c_k)R^k against the actual `tc_tm.TM` | PROVEN AS STATED |
| 3c | ι = 0; intersection with the natural enclosure | PROVEN AS STATED |
| 3d | mv3 float ι-centre | GAP (robustness, not triggered) |
| 4a | Merge validity checks at `c3beffc2` | GAP (checks the scope the driver declares, not the cover actually run) |
| 4b | A layer part never stands in for a missing or failed S4-3 part | PROVEN AS STATED |
| 4c | Provenance: lazy hashing; driver not pinned | GAP (not triggered) |
| 5 | Spot checks and recomputation | PROVEN AS STATED |

**Endgame factor structure** (exact, fmpq_mpoly through `tc_energy.Atoms`):
- SD = SΔ̃ = Dp·SDr;
- Bc = Dp²·Bcr;
- D2 = Dp³·D2r;
- D2r = c·dn·(e/Dp), with c free of κ, dn = (κ − A)(1 − V²), and κ-roots A
  and −F_N, both real;
- S, N and W divide D2r;
- κ-degrees: Dp 0, SDr 0, D2r 2, Bcr 3, dT 5.

So dT ≠ 0 at a strip κ gives every nonvanishing the endgame uses there:
- D2 ≠ 0 (closing the (χ, η) pair);
- b = Bc/D2 ≠ 0;
- dn, e ≠ 0 (ker restricted to Σ is span g);
- SΔ̃ ≠ 0 (P invertible);
- B0 ≠ 0 (commutation).

What remains depends only on κ, not on the piece (κ ≠ κ̄, Re σ(κ) < 0), or is
κ-independent (the gauge identity, and P4-3's 𝔴 > 0). S4-2 has its own pole
condition (RBAR = 2, tc_energy.py:52), and Theorem B′'s winding covers R.

**2f, the flaw.** The evaluation ball can miss the strip endpoints.
- `TM.eval` (tc_tm.py:118, 123–124, 139) and `LC.eval` (:270) evaluate on
  `arb((lo+hi)/2, (hi-lo)/2)`, and so does S4-3's `Bnd.eval_up`
  (tc_commutator.py:404).
- Because the midpoint is rounded to a float, this ball can miss lo or hi by
  about one ulp.
- On 200,000 random intervals, 55,581 missed. Dyadic bisection strips of
  [0, h] never missed (229,369 tested). The engine does pass non-dyadic
  strips, from `box_geometry` and the geometric Laurent splits.
- In practice Arb's radius round-up of up to 2⁻³⁰ relative covers the
  midpoint error, except on strips narrower than about 1e-7 relative. That is
  why no miss was observed.
- Fix: keep the float ball when it provably contains both endpoints, else use
  `arb(lo).union(arb(hi))`, and likewise for the recentred delta. Results stay
  bit-identical wherever no miss occurred.

**3d.** mv3's ι-centre ι_c = (1/t_lo + 1/t_hi)/2 is computed in floats
(tc_layer_mv.py:230). It falls outside the box's ι-ball for degenerate and
near-degenerate τ-ranges: 15% of ranges at relative width 1e-15, 0.18% at
1e-13, none at ≥ 1e-10. The covers never go below about 1e-9 relative, but
`piece_scale` does. Fix: take ι_c from the ball's midpoint.

**4a.** `layer_state` checks the scope fields. The driver writes those fields
from constants and the plan (s44_layer.py:278–284), not from the cover
actually run. The cover's own `la`, `lb`, `eps_cells` and `s_lo` are not
checked, and the driver hash (`driver_sha256`) is recorded but never checked.

**4c.** `loaded_code()` hashes sources at a worker's first part write, not at
import (s43_parts.py:197–210). An edit between import and first write would
be recorded as the code that ran. This did not happen in this run.

## 2. Conditions for the opt-in substitution rule

A passing `layer_t14` part d may carry plan piece p for every τ ≥ b1 only if
all of the following hold. Rbar and the real-pole bounds are not conditions.

1. **Eligibility.**
   - p.kind is tube or sonic, never centre.
   - p's S4-3 state is missing or stale.
   - An S4-3 part that failed or errored, paired with a passing layer part,
     goes to manual review, not to automatic substitution.
2. **Identity.**
   - The schema is `s44-layer/1`, d.id = p.id, and d.spec = spec_of(p).
   - For tube pieces: the record hash, the plan hash and
     d.inputs.cache_content_sha256 are equal, and d.inputs.cache = p.cache.
3. **Code.**
   - `code.files` is non-empty and every hash is current. No
     `--accept-code` for substituted pieces.
   - The required modules are present: tc_layer, tc_commutator, tc_tm,
     tc_energy and tc_reduce; lintube for tube pieces; tc_layer_mv if and only
     if `engine.mv` is set.
   - `git_dirty_src` is false.
   - `driver_sha256` equals the current hash of s44_layer.py, unless every
     check of item 5 is enforced.
4. **Result.**
   - status pass, `ok` true and `fails` empty.
   - certified_bound = result.certified_bound = result.worst, finite and
     < 0.
5. **The cover actually run**, read from `result.params`:
   - τ_min ≤ 14 and equal to scope.tau[0]; scope.tau[1] = "inf".
   - `eps_cells` contiguous from 0.0 to 6.0.
   - la = 0 and lb = p.h for tube pieces; la = p.lo and lb = p.hi for sonic
     pieces; scope.ξ consistent with these.
   - s_lo ≤ ξ(la)·τ_min, and 0 on the closing piece.
   - `laurent` equals (kind == sonic).
   - The evaluator is sigma_box, sigma_box_mv3 or sigma_box_mv, with mv only
     on tube pieces.
   - `combine` equals `engine.combine`, and is false for sonic pieces.
6. **Global.**
   - checks.json is valid for the current code. E1, E3 and E5 justify the
     sonic raises; E3 also justifies the (tt·dT_k)′ raise.
   - b1 = 14, c = −1/32 and EPS_MAX = 6 are unchanged in tc_commutator and
     tc_layer.
7. **Assembly.**
   - The substituted interval enters `direct_iv` and `both_iv` and leaves
     `missing`; its code hashes join `code_union`.
   - `complete` and `strip_exclusion_certified` are recomputed.
   - The detail entry records that layer_t14 carries both τ ranges, and that
     Bcr, D2r ≠ 0 at strip κ follow from dT ≠ 0, in place of Rbar and the
     real-pole bounds.

## 3. Computations

1. **Tests:** `tests/test_tc_layer.py`, `test_s43_merge.py` and
   `test_s44_driver.py` gave 58 passed and 4 skipped (the P4_SLOW set) in
   24 s.
2. **Exact factorizations** (fmpq_mpoly, about 10 s each): as in §1.
3. **Identity checks** (sympy):
   - the Σ residual, the B1 split, and the ν and ρ normalizations are all 0;
   - two perturbed controls are nonzero.
4. **Arb semantics:** division by a zero-containing ball gives nan, 0/z
   included, and nan·0 = nan.
5. **`TM.eval` ball test:** 229,369 dyadic strips gave 0 misses; 200,000
   random intervals gave 55,581.
6. **ι-centre test** at precision 320: the miss rates of 3d.
7. **Independent point evaluation, 231 s.**
   - Path: component values → constant TMs → `AtomKP` → `build_pairs` → the
     reviewer's own S chain.
   - It reproduces the note's §6 values −3.2098, −15.457, −1.141 and −42.956.
   - 16,059 points on 7 pieces: grids with τ up to 1e10, random points, and
     coordinate ascent for the maximum. No midpoint lies above any bound.

   | Piece | Float sup | Bound τ ≥ 14 (ratio) | Bound τ ≥ 224 (ratio) |
   |---|---|---|---|
   | deep/deep-128 | −0.043285 | −0.021657 (0.500) | −0.021728 (0.502) |
   | tube/full-244 | −0.418746 | −0.209392 (0.500) | −0.210185 (0.502) |
   | gap/gap-000 | −1.2061 | −0.603886 (0.501) | −3.89314 (0.565 of −6.893) |
   | sonic/s40 | −0.938416 | −0.469208 (0.500) | −0.757218 (0.664 of −1.140) |
   | sonic/s00 | −0.8978 | −0.500341 (0.557) | −0.763316 (0.850) |
   | tube/full-138 | −6.1726 | −3.37892 (0.547) | −5.93162 (0.949 of −6.247) |
   | tube/full-181 (fresh) | −2.9444 | −1.47316 (0.500) | −1.86487 (0.623) |

8. **sonic/s00 near ξ = 0:** 1,600 points at ξ ∈ [1e-7, 1e-6]. The supremum
   is −0.8978; no tight point lies above −0.5003 (τ ≥ 14) or −0.7633
   (τ ≥ 224).
9. **Recomputation from scratch** with `TM.eval` instrumented. Every count
   matches the stored part exactly:

   | Piece | τ ≥ 224 bound | τ ≥ 14 bound | Instrumented strips, misses | Time |
   |---|---|---|---|---|
   | tube/full-244 | −0.21018543312934612 | −0.20939233900662965 | 85,680, 0 | 262 s |
   | sonic/s40 | −0.7572182020906124 | −0.4692083392578729 | 229,222, 0 | 113 s |
   | sonic/s00 | −0.7633164782585093 | −0.5003413981238266 | 10,878, 0 | 82 s |
   | tube/full-181 | −1.864870243105424 | −1.4731633357498715 | 10,416, 0 | 98 s |

   Point values against both enclosures on 1,380 test boxes: none disjoint.
10. **Part scan:** all 451 stored layer parts satisfy the conditions of §2,
    item 5.
11. **Modification times:** no `src/p4` file changed after the workers
    started.

## 4. Wording fixes, must-fix items, ledger level

**Wording fixes.**
- `s4-4-layer-engine.md` §1 and §13.2: replace the Rbar and real-pole claims
  with the factor argument of §1 above.
- §4: every ratio has k + mN − mD = 0, except λ²·nV′/dT, which carries ℓ^{+1}.
- §5 (v): add the evaluation-ball caveat until it is fixed.
- §10: the ξ and ε scope fields are declared by the driver.
- §12: state the TM semantics: valid at each point of [0, R] with R = h, no
  separate remainder term, eps_z in c₀.
- `s4-commutator.md` §3 (:125–129) and §10 (:334–339): qualify "certified per
  piece" for the substituted pieces.

**Must-fix before the layer results enter the ledger.**
1. Implement the opt-in rule with every condition of §2, with a test for each
   new rejection mode.
2. Correct the written justification as above.
3. Fix `TM.eval` and `LC.eval`, and either recompute the layer parts or accept
   the new code citing §3. Flag the same flaw to the S4-3 owner
   (`Bnd.eval_up`).

**Should-fix.**
- ι-centre from the ball (tc_layer_mv.py:230).
- Hash sources at import and re-verify at write (s43_parts.py:197–210).
- Pin s44_layer.py in the merge.

**Ledger level.** With the must-fix items done, the strip exclusion via the
layer engine merits CERTIFIED.
- The decisive step, S < 0 and dT, dB, nB ≠ 0 on every piece × strip, is
  machine-checked.
- The added lemmas are routine: exact algebra, the (tt·dT_k)′ analyticity
  raise, and the mean value theorem.
- The row inherits S4-3's prose lemmas: the energy identity, the sonic
  boundary terms and the endgame. S4-5's review must accept them as routine,
  or both paths drop to WRITTEN_PROOF.
- The strip exclusion is not complete until the 30 centre pieces land.

## 5. Coordinator's note (2026-10-17, after the review)

The float-midpoint ball of 2f is not confined to the files the review read. A
repository-wide search for `arb((a+b)/2, (b-a)/2)` found it only in P4:

| File | Occurrences | Used by |
|---|---|---|
| `tc_commutator.py` | 10 | S4-3 |
| `tc_tm.py` | 3 (plus the recentred delta) | S4-3, S4-4 |
| `tc_energy.py` | 1 | S4-2 (P4-8) |
| `tc_monotone.py` | 1 | P4-3 |
| `modecount.py` | 1 | Theorem B, P4-6, an input of Theorem C |
| `scripts/s43_centre_tq_diag.py` | 1 | diagnostic only |

The decision is to fix every occurrence the same way, keeping the float ball
when it provably contains both endpoints so that unaffected results stay
bit-identical. Then re-run every affected certificate on the fixed code and
compare against the stored outputs:
- the layer parts;
- checks, Lemma D and P4-3's centre;
- Theorem B's winding, analyticity and Krawczyk;
- `tc_energy`;
- `tc_monotone`.

Earlier timings put the non-layer re-runs at about 1.7 core-h. None of these
certificates records code hashes, so the comparison is by output.
