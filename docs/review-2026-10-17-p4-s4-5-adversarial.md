# Independent adversarial review: P4 S4-5, the Theorem C chain (all but the centre computation)

Date: 2026-10-17. Reviewer: an independent agent with no part in S4-1 … S4-5.
Reviewed: branch `p4/theorem-c-s4-3` at `75869b91` (S4-5 assembly note
`problems/P4/notes/s4-5-assembly.md`, sha256 `99eef2ba…`; HEAD `ee7a2625` adds
only a timing bench record, not read by the merge).  `P4/` means
`problems/P4/`.

How the review ran:
- Read-only in the worktree.  Every computation ran in a `git archive 75869b91`
  scratch copy (the gitignored tube caches copied in), at nice 10, at most two
  of my processes at once next to the two benchmark processes, about 0.15 core-h
  in all.
- The cover and provenance items were reviewed after `a0a1b5ff` and `ccc74f3a`
  landed and the coordinator declared the files stable.
- My scripts and their outputs (list in §9) are committed with this document
  in `problems/P4/results/theorem_c/review_s45/`.  Several verdicts below say
  "closed by this review's check"; those checks are the evidence, because the
  notes assert the steps they back without proof.

**Overall verdict.**
- **No FLAW.** I found nothing wrong in the chain.  The S4-3/S4-4 argument and
  S4-2 (P4-8) are sound for every κ they claim.  The centre computation is out of
  scope here, and Theorem C is still missing it.
- **Five GAPs**, all closed or closable at negligible cost; two of them were
  asserted in the notes from float evidence or a mis-cited certificate:
  1. *S8, sonic exponent.*  The residues behind γ > 1 are cited from the S4-2
     closing box, which is valid only for Re κ ≥ 6, and from large-κ quotient
     hats, which are not fixed-κ residues.  Closed here: the residue identities
     hold as κ-polynomial identities, in balls, at the sonic point, and an
     integer-rigidity lemma makes them exact.  γ ≥ 1.0616718.
  2. *S8 and S4-2, centre limits.*  T∞ = −3, B0 → 0 and tt → 0 are asserted
     from an exponent heuristic and floats.  E2 gives 𝒬 → 9/4, which fixes only
     the exponent difference, not T∞.  Closed here by exact centre-layer
     identities: T = −3 + O(t²), B0 = O(t²), tt = O(t²) for every κ ∉ {0, 1}.
  3. *G1, the E-form outside R.*  The note understates it.  Besides the sonic
     order-0 solve, the centre regular family (certificate (C) of
     `s2-theorem-b.md` §3.5) is also certified on R only.  Both halves are
     settled here: det M₀(κ) = d₁(κ − κ\*) with κ\* = −0.0990290820550 ± 3e−14,
     and det(nP₀ − G₀(κ)) is κ-free by an exact zero pattern.
  4. *The S4-2 x-cover.*  The join at x = −4.5 between the tc_energy tube region
     and Lemma D is recorded nowhere.  It holds with an exact overlap of only
     4.6e−17.
  5. *The class wording.*  P4-6's class carries the normalization A_p(0) = 1
     and a radius; the energy arguments use the larger, un-normalized class.
     The two are equivalent on Re κ > κ\*, by item 3.
- **Ledger.**
  - P4-8: "conditional on Lemma D" can be dropped, so the row is CERTIFIED, once
    the checks of items 2 and 4 are committed.
  - P4-8's "hence E(κ) ≠ 0" clause holds as written only with the G1 lemma (§2).
    Until then, restrict it to [6, 15] × [−14, 14].
  - P4-9 and P4-10 stay CONJECTURED until the 30 centre pieces land.  I accept
    S5–S8, S12 and S13 as routine, given the committed checks, so both rows then
    merit CERTIFIED.
- **Cover and provenance** (after `a0a1b5ff`, `ccc74f3a`).
  - The merge's exact cover is correct.  I reproduced it with independent code.
  - The two provenance fixes work end to end.  Two parts I recomputed in the
    archive copy were bit-identical but rejected by the strict rule for missing
    git provenance.

## 1. Verdicts

Basis: **X** exact identity checked by machine; **M** ball/interval
certificate; **W** written argument.  "Level" is what the item supports after
this review.  "+commit" means the review's check must be committed for that
level.

| # | Item | Verdict | Basis | Level |
|---|---|---|---|---|
| 1a | Class identical across P4-6, P4-8, S4-3, S4-4 | GAP (wording): P4-6 normalizes A_p(0) = 1 with radius ≥ 0.065; S4-2/S4-3 use nonzero p with N_p(0) = 0. Equivalent on Re κ > κ\* (2c) | W + M | CERTIFIED-grade (+commit `g1_order0`) |
| 1b | Eigenvalue form of Theorem C + B′ on R give the mode count without E outside R | PROVEN AS STATED | W | routine |
| 1c | G1: E-form on Ω | GAP, closable (§2) | M + X + W | E-form: CERTIFIED-grade after the §2 lemma and +commit |
| 1d | Higher-order sonic resonances in Ω | PROVEN AS STATED: det(nD + E(κ)) = det D · n³(n − σ(κ)), E exactly rank one, Re σ ≤ −0.0901 on Re κ ≥ 0 | X + M | CERTIFIED-grade |
| 1e | Each energy argument applies to every solution of the class, whatever the Frobenius structure | PROVEN AS STATED (§2d) | W | routine |
| 2 | Region algebra Ω ⊂ {Re κ ≥ 6} ∪ strip ∪ conj(strip), closed sets as certified | PROVEN AS STATED | W | routine |
| 3a | S6, the χ₁ energy identity | PROVEN AS STATED (sympy, residual 0; negative controls nonzero) | X | CERTIFIED-grade (+commit `check_s45_identities.py`) |
| 3b | `direct_S_boxes` and the layer Σ evaluate the same S | PROVEN AS STATED (code read; Σ = Sλ²/τ² by sympy) | X + W | CERTIFIED-grade |
| 3c | Elimination, normal form, commutation (S4, S5, s4-energy §1) and the polynomial forms TN1, U1, QNm | PROVEN AS STATED (sympy) | X | CERTIFIED-grade |
| 4a | S7, the one-sided limit argument | PROVEN AS STATED | W | routine |
| 4b | S8 sonic: γ > 1 for ε ∈ [0, 6], c = −1/32, from certified quantities | GAP (citation) → closed here | M + W | CERTIFIED-grade (+commit `sonic_residues`) |
| 4c | S8 centre exponents | GAP → closed here (exact) | X | CERTIFIED-grade (+commit `centre_layers`) |
| 4d | S4-2 boundary terms (`s4-energy.md` §2) | sonic: PROVEN AS STATED for Re κ ≥ 6; centre: GAP (T∞ = −3) → closed by 4c | M + X + W | CERTIFIED-grade (+commit) |
| 5a | Every polynomial needed nonzero at a strip κ divides dT up to Dp powers | PROVEN AS STATED (exact; D2r = 16S⁴W²·dn·(e/Dp)) | X | CERTIFIED-grade |
| 5b | Import-time factorizations `tc_energy.py:79–121` | PROVEN AS STATED (sound exact divisions; re-run) | X | CERTIFIED-grade |
| 5c | dT_top ≠ 0 and dB_top ≠ 0 enforced | PROVEN AS STATED, per piece (stronger than per box) | M | CERTIFIED-grade |
| 5d | S12, the endgame | PROVEN AS STATED | W + X | routine |
| 6 | S13, conjugation | PROVEN AS STATED (notation clash, §6) | W | routine |
| 7a | S4-2 chain (a)–(d) ⟹ no eigenvalue | PROVEN AS STATED | W | routine |
| 7b | Lemma D from `parts/lemma_d_certified.json` | PROVEN AS STATED (24 recorded hashes current; E2 valid; y → 0 included; the mC margin is conservative) | M + X | CERTIFIED |
| 7c | S4-2 x-cover, including the join at −4.5 | GAP (unrecorded) → checked here: overlap 4.6e−17 | X | CERTIFIED-grade (+commit `s42_join.py`) |
| 7d | Drop P4-8's "conditional on Lemma D" | YES, with 4c and 7c committed | — | CERTIFIED |
| 7e | P4-8's "hence E(κ) ≠ 0 on {Re κ ≥ 6}" | GAP (needs G1) | — | restrict to [6, 15] × [−14, 14], or add the §2 lemma |
| 8a | Exact contiguity of cache steps | PROVEN AS STATED (merge and independent code) | X | CERTIFIED |
| 8b | Cross-cache and sonic joins | PROVEN AS STATED | X | CERTIFIED |
| 8c | The overlap the centre must meet | PROVEN AS STATED: 0.06 in x (c29 reaches −4.5 − 3.7e−17; the deep end is −4.56) | X | — |
| 8d | Provenance fixes `a0a1b5ff` | PROVEN AS STATED (tests; end-to-end rejection) | M | — |
| 8e | Merge output wording | GAP (cosmetic): `xi_open_max = Infinity`, `full_strip_tau = [14, 224]` remain | — | — |
| 9 | Levels proposed in assembly §6 | Agreed, with the conditions of §11 | — | — |

## 2. G1: the E-form outside R

**(a) What the eigenvalue form needs.**  S4-2 and S4-3 exclude any nonzero
solution p that is
- analytic at x = 0;
- on Σ;
- regular at the centre (analytic in t, scaled variables bounded);
- subject to N_p(0) = 0.

They never use E, p̃, r₁ or r₂.  Theorem B′ gives "E = 0 ⟺ eigenvalue" on R.
So:
- the eigenvalues with Re κ ≥ 0 are those in R;
- these are κ̄ and κ₁;
- the mode count follows.

E outside R is not needed (assembly §4.3 A4 is right).

**(b) The sonic half.**  `g1_order0.py` builds `_run4`'s order-0 matrix exactly
as `linsonic4.py:142–150` does, on the certified V₀ box (A1 background,
256 bits).  det M₀ is affine in κ, since only the ℓ·G row depends on κ:

- det M₀(κ) = d₀ + d₁κ, with
  - d₀ = 1.17284764762 ± 4.4e−12;
  - d₁ = 11.8434668208 ± 2.8e−11.
- Its root is **κ\* = −0.0990290820550 ± 2.9e−14** (real).

So det M₀ ≠ 0 on Re κ > κ\*, in particular on all of Re κ ≥ 0.

The note's "singular iff σ(κ) = 0" is also confirmed:
- σ(κ) = σ_a + σ_bκ, with σ_a = −0.0901059705079 ± 5e−14 and
  σ_b = −0.909894029492 ± 6e−14;
- its root equals κ\* (d₀σ_b − d₁σ_a ∈ [± 1.9e−13]).

The structural reason is that the constraint gradient annihilates the right
null vector of the fluid block P_fl(u₀): C̃·r_null ∈ [± 1.0e−14].
Then det M₀ = A₀·ℓG₀(κ)·(0, 0, C̃_V, −C̃_W) is proportional to
σ(κ) = ℓG₀(κ)·D⁻¹e₄.

**(c) Higher orders, and the class equivalence.**
- Rows 0–2 of E₀ and E₁ are exactly zero, so E(κ) has rank ≤ 1, and det D =
  −9.60210331438 ≠ 0.  Hence
  det(nD + E(κ)) = det D · n³ (n − σ(κ)) identically.
- Re σ(κ) ≤ σ_a < 0 on Re κ ≥ 0, so M_n is invertible for every n ≥ 1.  No
  resonance occurs anywhere in Ω, real κ > 15 included.
  - There σ(κ) is a negative integer (σ = −14 near κ ≈ 15.29).
  - That affects only the x^σ solution, which is unbounded.
- The admissible family {analytic at 0, c(0) = 0, N_p(0) = 0} is
  ker(ℓG₀(κ)) ∩ {N_p = 0} ∩ {c(0) = 0}, which is one-dimensional.  Its
  A_p-component is ∝ det M₀ ≠ 0.
- Hence, for every Re κ > κ\*:
  - the normalized class of P4-6 and the un-normalized class of S4-2/S4-3
    coincide (item 1a);
  - p̃ exists and is unique.
- Convergence of the formal series is classical (regular singular point of the
  first kind: p′ = x⁻¹B(x)p with B = x·adj P·G/det P analytic).
- The radius 0.065 of P4-6's class is a property, not a restriction:
  analyticity at 0 is enough.

*Failure scenario of the class wording (1a).*  Suppose det M₀ had vanished at
some κ ∈ Ω.
- P4-6's normalized class would then be empty there, while a solution with
  A_p(0) = 0 could still exist.
- A row phrased "in the class of P4-6" would assert nothing at that κ, while
  the energy arguments prove the un-normalized statement.

*Fix:* define the class without normalization, as in §11.  By (b) the two
classes agree anyway on Re κ > κ\*.

**(d) Any Frobenius structure is harmless.**
- Since Re σ < 0 on Re κ ≥ 0, every solution that is bounded (a fortiori C¹)
  as x → 0⁻ is analytic there.  The other branch is ~ |x|^{Re σ} → ∞.
- The energy arguments use only:
  - (i) analyticity of p at 0 (χ = O(x), and χ₁, χ₁′ bounded);
  - (ii) analyticity in t at the centre;
  - (iii) the coefficient residues and limits of §4.
- Items (i) and (ii) hold for every solution in the class.  Item (iii) is a
  property of the equation, not of the solution.

**(e) The centre half.**  The E = 0 ⟺ eigenvalue equivalence of
`s2-theorem-b.md` §3.8 also uses (C): the exponent-0 family is exactly
span(r₁, r₂).  That too is certified only on R's boxes.

`centre_pattern.py` shows exactly, as polynomials in (n, w, v) at t = 0 with
θz = 0 (`linscaled.full_system`):
- P₀ has rows 1 and 3 equal to (0, 1, 0, 0) and (0, 0, n, 0);
- G₀ has rows 1 and 3 identically zero;
- the κ-part of G is the single entry G[2][2] = −3κ.

Expanding along rows 3 and 1 removes column 2, so
det(nP₀ − G₀(κ)) = ±n · (n∞n) · (n + 3) · (P₂₃n − G₂₃) is κ-free, with
P₂₃ = 28.79 and G₂₃ = −86.38.  The ball values
(`g1_centre.py`, certified (a\*, μ\*) series) are identical at κ ∈ {0, 1, 7,
3 + 50i, 1000i}: −569.685, −3560.53, −11536.1, … for n = 1 … 6.  The exponents
are {0, 0, −3, −3} for every κ, and G₀(κ) has rank 2 with the κ-free
(Â, ñ)-minor diag(−3, 28.24).

The remaining ingredients of §3.8 are already κ-independent:
- Σ invariance (exact);
- the constraint exponents (ball);
- e_W ∉ Σ(x_d);
- independence of r₁(x_d), r₂(x_d).

**Failure scenario the gap allowed.**  Two things could have broken the E-form:
- if κ\* had been a real number > 15, E (normalized by A_p(0) = 1) would be
  undefined there;
- if nP₀ − G₀(κ) were singular for some n ≥ 1 and κ ∈ Ω, the regular family
  could be 3-dimensional and the equivalence would fail.

Then "E(κ) ≠ 0 on Ω" would be meaningless or false, even though no eigenvalue
exists.  Neither happens.

**Recommendation.**
- State P4-9, P4-10 and P4-8's main clause in **eigenvalue form** now.  That is
  all the mode count needs.
- **Close G1 as a corollary.** It is cheap:
  - commit `g1_order0.py`, `g1_centre.py`, `centre_pattern.py` and their outputs;
  - write a half-page lemma combining (b), (c) and (e) with classical Frobenius
    convergence and the κ-independent ingredients of §3.8.
- After that, the E-form "E(κ) ≠ 0 on Ω" holds as FORMULATION §2 states it, and
  P4-8's E-clause can stay.  Until then, restrict P4-8's E-clause to
  [6, 15] × [−14, 14], where it follows from B′.

## 3. Region algebra

- Ω = {Re κ ≥ 0} \ R with R = [0, 15] × [−14, 14] closed.
- A point of Ω has either Re κ > 15, or 0 ≤ Re κ ≤ 15 with |Im κ| > 14.
  - Re κ ≥ 6 lies in S4-2's closed half-plane.  B₀ and B₁ are suprema over
    Re κ ≥ c₀, by maximum modulus on the line, and ε = c₀ is included.
  - Re κ ≤ 6 with Im κ > 14 lies in the strip [0, 6] × [14, ∞).  The ε-cells
    tile [0, 6] exactly (checked in all 546 parts), and τ = 14 and τ → ∞
    (pseudo-boxes) are included.
  - Im κ < −14 lies in the conjugate strip.
- Every boundary is covered by a closed certified set.  Nothing needs an open
  boundary.
- The apparent singularities sit at real κ and do not touch Ω:
  - κ = A(x) lies in (1, 1.89];
  - κ = −F_N(x) lies in [0.3549, 1) (float, `s4-reduction.md` §1.4).  The
    note's "[κ̄, 1.89]" misses the dip of −F_N below κ̄.  This is cosmetic.
- σ = 0 at κ\* = −0.0990 < 0 (§2b).

## 4. The strip identity, the code paths, and the boundary terms

**S6.**
- The hand derivation in assembly §4.1 is right.
- `check_s45_identities.py` (sympy, no repository code) checks I5, a pointwise
  differential identity:
  BT′ = g|u′|² + g ν Im(ūu′) − g K |u|²
  for every solution of u″ + A1u′ + B1u = 0, with
  - g′ = ρg, ρ = Re A1 − c·tt, ν = Im A1 + τ·tt;
  - f = (c + τ²/c)tt − (τ/c)ν;
  - K = Re B1 − (τ/c)Im B1 − (f′ + ρf)/2;
  - BT = Re[(1 + iτ/c) g ū u′] + g f |u|²/2.
- The residual is 0.  Changing the sign of the (τ/c)Im B1 term gives a nonzero
  residual (control).
- The same script checks:
  - I1, the elimination χ″ − Tχ′ − Uχ = 0 with T = a + d + b′/b and
    U = a′ + bc − a(T − a);
  - I2, the Liouville normal form;
  - I3, the polynomial forms of `tc_energy.py:100–108`:
    - T = TN1/(SD·Bc·D2);
    - U = U1/(SD·Bc·D2²);
    - 𝒬 = QNm/(SD³Bc²D2²), with D1 = SΔ̃·d/dx;
  - I4, the commutation as an operator identity;
  - I6, the layer normalization of `tc_layer.sigma_box`, transcribed term by
    term: Σ − Sλ²/τ² = 0 (control with one sign flipped: nonzero).
- I3c was also spot-checked numerically with concrete functions; the sides
  agree to 1e−15.

**The code paths.**
- `direct_S_boxes` (`tc_commutator.py:1314–1507`) evaluates, box by box, the
  definitions A0 = nT/dT, A1 = A0 − (log B0)′, B1 = B0 + A0′ − A0(log B0)′,
  ν, ν′, ρ, f, f′, K and S = K + ν²/4.
  - The w-Horner hats return P(κ)/κ^deg.
  - Their ratios restore nT/dT and rU/dB exactly.
- The derivative slots are x-derivatives in every region: X′ = D1(X)/SΔ̃
  (`tc_tm.AtomKP.jet_kp`), whatever the TM variable.
- `tc_layer` builds V, Y, (log B0)′, (log B0)″ and W from the same
  `build_pairs` polynomials.  Its Σ is Sλ²/τ² by I6.  The `combine` and mv3
  variants rewrite the same quantities.
- `point_S_direct` is the same chain as `direct_S_boxes`.
- So both engines bound the same S.

**S7.**  This is correct, and needs only BT → 0 at both ends.
- ∫_{x₁}^{x₂} g|S||u|² ≤ BT(x₂) − BT(x₁) by Young.
- The left side is monotone in the interval.
- g > 0, and S ≤ bτ²/ξ² < 0 pointwise.

**S8, sonic end (GAP 1, closed).**  The γ formula needs three fixed-κ
quantities at x = 0:
- R_T(κ) := lim xT;
- R_B(κ) := lim x²B0;
- Res tt.

The assembly cites two sources, and neither gives them on the strip:
- the S4-2 closing box `xT = τ₀ + τ₁κ ± 2.3e−5`.  `closing_box` bounds the
  remainder by `line_sup` over **Re κ ≥ c₀ = 6** (`tc_energy.py:704`).  The
  stored values are float midpoints plus a 1.75e−5 bound, not 2.3e−5.
- the S4-3 "u-hats".  These are ξ → 0 hats of the **large-κ** quotient
  coefficients.  The note's own §6 (the b0 ξ⁻³ hat cancelled by the remainder
  at fixed κ) shows that those limits need not commute with the fixed-κ
  residue.

*Failure scenario.*  Suppose the fixed-κ residue of A0 differed from
−(1 + σ(κ)) by about 0.06 somewhere on the strip, for instance through a
remainder term invisible to the large-κ hats.  Then γ ≤ 1 at ε = 0, and
BT = O(|x|^{γ−1}) need not vanish.  The cited data cannot exclude this.

*Closure.*  `sonic_residues.py` evaluates, at the sonic point (balls, V₀ box,
E_D = D′(0) = 7.56057255908 ± 5e−12), the κ-polynomials of
- TN1r(u₀; κ) − (1 + σ(κ)) E_D SDr D2r Bcr(u₀; κ): all seven coefficients
  contain 0 (radii ≤ 3.5e−9, relative ≤ 1.7e−12);
- −(dU6/dx)(0; κ) − (1 + σ(κ)) E_D³ SDr D2r² Bcr(u₀; κ): all coefficients
  contain 0.  U6(u₀) ∋ 0 coefficientwise, as E5 says.

It also gives Res tt = TN1r₆/(E_D SDr D2r₂ Bcr₃) = −0.909894029492 ± 6e−14,
which is −(1 − σ₀).

Exactness follows by integer rigidity:
- For κ with σ ∉ ℤ, the scalar χ-equation is Fuchsian at 0.  Its solution
  space is spanned by
  - χ_a = ℓ(p_a), p_a an analytic Σ-solution not ∝ g, with integer exponent
    m_a ≥ 1, because Δ̃(0) = 0;
  - χ_s = ℓ(p_s), p_s = x^σφ, with exponent 1 + σ + j, j ≥ 0;
  - each is ≢ 0: otherwise η ≡ 0 too (b ≠ 0), and the endgame argument puts
    the solution in span g, which is analytic.
- So R_T − σ = m_a + j is an integer on an open set.  It is a rational function
  of κ, hence constant.  The ball check pins it to 1.
- Then m_a = 1, j = 0 and R_U = −(1 + σ).  So R_T = R_B = 1 + σ(κ) exactly, and
  Res tt = −(1 − σ₀).

This gives:
- Res A1 = 2 − R_T = 1 − σ(κ);
- **γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.061671782 ± 9e−11** for ε ≥ 0;
- the sample balls agree, e.g. γ(3 + 20i) = 3.7913539.

ν is bounded at 0, because the residues of Im A1 and τ·tt cancel; f = O(1/|x|).
So BT = O(|x|^{γ−1}) → 0.

**S8, centre end (GAP 2, closed).**
- The assembly derives T → −3 from "χ-exponents {0, −3}" and says E2 backs it.
  E2 gives only 𝒬 → 9/4.  In the normal form, 𝒬∞ = (r₁ − r₂)²/4 fixes only
  the exponent difference, not T∞ = r₁ + r₂.
- *Failure scenarios.*
  - The lowest centre layer of dT could vanish at some κ ∈ Ω.  Then A0 has a
    pole in t, g ~ exp(±C/t), and BT need not vanish.
  - In S4-2, T∞ = +3 would give f ~ e^{3x/2} → 0, and v = χ/f would grow.
- *Closure.*  `centre_layers.py` works exactly in fmpq_mpoly, in E2's grading
  p = −e_N + 2e_W + e_V with A → 1.  It finds:
  - val(TN1r) = val(denT) = 7 and **val(TN1r + 3·denT) = 9**;
  - val(denU) = 8 and **val(U6) = 10**;
  - val(TN1r₆) = 9 > val(denT₅) = 7.
- The lowest layers are
  - denT: −12288 n⁵w⁶ κ(κ−1)⁴;
  - denU: −786432 n¹⁰w⁹ κ(κ−1)⁶.

  Both are nonzero at the centre for κ ∉ {0, 1}.
- Hence, for every κ ∈ Ω: T = −3 + O(t²), A0 → 3, B0 = O(t²), tt = O(t²).
- `centre_U6_roots.py` adds a float diagnostic.  The t¹⁰ layer of U6, at the
  certified centre values (2nv + 1 ∋ 0), has κ-roots 4.023, −5.877, −0.598,
  1.119 and five more within 0.002 of κ = 1.  None is in the strip or in
  Re κ ≥ 6.  So B0 = βt²(1 + O(t)) with β ≠ 0 there: k = 2.
- The regular χ = c₀ + O(t²) (the t¹ coefficient solves 4c₁ = 0).  More
  generally χ₁, χ₁′ = O(t^k) and g ~ t^{3−k}, so BT = O(t^{3+k}) → 0.  This
  holds for any k ≥ 1, so the float value k = 2 is not needed (it gives
  O(t⁵)).

**S4-2's boundary terms (`s4-energy.md` §2).**
- *Sonic end:* correct for Re κ ≥ 6.  The closing-box enclosure applies there,
  and §4's rigidity lemma makes it exact.  |v̄v′| = O(|x|^{−Re σ}) and
  −Re σ ≥ 5.55.
- *Centre end:* needs T∞ = −3 (now exact).  Then f ~ e^{−3x/2}, v, v′ =
  O(e^{3x/2}), and [v̄v′] → 0.

## 5. Factorizations, nonvanishing and the endgame

- `factor_checks.py` (exact, through `tc_energy.Atoms`):
  - e = Dp·(e/Dp), with e/Dp = 4SW(Q₁ + κN) and N | Q₁, so
    Q₁ + κN = N(F_N + κ);
  - **D2r = 16S⁴W² · dn · (e/Dp)**;
  - SDr = 4S²W (κ-free), SΔ̃ = Dp·SDr = det P;
  - D2 = Dp³D2r and Bc = Dp²Bcr;
  - κ-degrees of (Dp, SDr, D2r, Bcr): (0, 0, 2, 3).
- Every quantity the strip endgame needs nonzero follows from dT ≠ 0 at a strip
  κ:
  - D2 (closing the pair);
  - b = Bc/D2;
  - dn and e (the kernel lemma);
  - g_N = N′ + κN;
  - SΔ̃ (P invertible).

  B0 ≠ 0 follows from nB, dB ≠ 0.
- **`tc_energy.py:79–121`** is sound.  `_dpow` divides exactly (fmpq_mpoly
  `divmod` remainder zero) and raises unless the multiplicity is the expected
  one.  U1 = Dp⁶U6 is checked as an identity.  The checks run at every
  `Atoms(full=True)` construction, and passed in my runs.
- **dT_top and dB_top** are inverted per piece by `TM.inv`
  (`tc_tm.py:152–183`).  It raises unless |q| < 1 on the whole TM domain, so
  the top is nonzero on the whole piece:
  - in the layer engine at `tc_layer.py:355` (tt) and `:379` (u2);
  - in the S4-3 path through `inv_of` in `analyze_piece`.
- The kernel lemma (dn·e ≠ 0 ⟹ ker(ℓ, m)|_Σ = span g) was re-derived: three
  independent covectors.  The endgame of assembly S12 is correct:
  - p = λg;
  - P(λ′g) = 0, so λ′ = 0;
  - λ N₀(κ − κ̄) = 0 with κ̄ = −F_N(0) = 2 − A₀ + 2W₀/3.

## 6. Conjugation (S13)

The lemma is correct and complete.
- P, DQ, Ψ and P_s have rational coefficients and the background is real, so
  G(κ)‾ = G(κ‾).
- Each class condition is preserved:
  - (i) and (iii): conjugate power series in x, resp. t, with real scalings;
  - (ii): c_κ(x, p)‾ = c_{κ‾}(x, p̄);
  - (iv): N_p̄(0) = 0.
- *Notation.* S13 (and `s2-theorem-b.md:535`) uses κ̄ for the complex
  conjugate, while κ̄ is also the gauge eigenvalue in S1, S12 and §1.2.  Write
  κ\* or κ^c for the conjugate.

## 7. S4-2 / P4-8

**Chain (a)–(d).**  The chain is sound.
- (b) and q₂ > 0 give Im𝒬/τ ≥ 2εq₂ + q₁ − B₁ > 0 for ε ≥ c₀.
- For τ = 0, φ(ε) = q₂ε² + q₁ε + q₀ is increasing on [c₀, ∞) because
  φ′ ≥ 2c₀q₂ + q₁ > B₁ ≥ 0.  So 𝒬(ε) ≥ φ(c₀) − B₀ > 0.
- R′ is analytic on Re κ > 2 (poles A, −F_N and the Bcr roots, all with
  Re ≤ RBAR = 2) and R′ → 0, so maximum modulus applies on the line and the
  Cauchy shift (DELC = 1.2) is legal.
- (d) also gives b ≠ 0 and dn·e ≠ 0, and P4-3 gives P invertible.

**Lemma D (`parts/lemma_d_certified.json`, sha256 `41fda13b…`).**
- `ok`, 1100/1100 subranges, the first [0, ·], so y → 0 is included.
- The certified (a\*, μ\*) balls are used.
- E2 = (14, True) holds, so val(Z∘u) ≥ 16 against val(QD) = 14.
- All 24 recorded source hashes equal the current files.
- The margins are hat-level:
  - (b): 2c₀q₂ + q₁ − B₁ = y²·mB;
  - (c): the code's mC = H + 9/4 with H the y²-hat.  The true condition is
    9/4 + y²H > 0, and it follows from mC > 0 because y² < 1, so the check is
    conservative.
- Hurwitz runs on the Bcr hat cubic; poleA and poleFN confirm A, −F_N < 2.
- `parts/replays/lemma_d_certified.json` records pre-ball-fix hashes for
  `tc_commutator`, `tc_energy` and `tc_tm`.  It is a historical replay; cite
  the certified file.

**The x-cover of S4-2 (GAP 4).**  No merge checks it.  `s42_join.py` checks it
exactly:
- tc_energy's tube region keeps the steps with float(x) − h ≥ −4.5 (274
  steps).
- Its lowest step ends at x = −4.500000000000000083266727 (exact dyadic,
  radius 0).
- Lemma D covers y ≤ YCUT = 0.011108996538242306, which is **below** e^{−4.5},
  i.e. x ≤ ln YCUT = −4.500000000000000037288925.
- The regions overlap by **4.6e−17**.  There is no gap, but it holds by
  rounding luck.
- The other joins are exact:
  - the tube's top step x₀ is exactly the double −0.05, the start of the sonic
    ladder;
  - the ladder ends at the double −1e−8, the start of the closing box.
- The recert run used the same cache (`tube_info` matches).

*Failure scenario.*  A tube rebuild whose step sum rounds the other way, or a
YCUT rounded down further, would open a sub-1e−16 sliver on which (a)–(d) are
unchecked, silently.

*Fix.*  Commit the check.  Better, make the regions overlap by a margin: run
Lemma D to a ycut slightly above e^{−4.5}, or keep one more tube step.

**Answers.**
- "Conditional on Lemma D" can be dropped, given 4c and 7c (committed).
- The E-clause: see §2.
- Evidence must cite `parts/lemma_d_certified.json` and the tube, sonic and
  closing blocks of `tc_energy.json` (with `recert_2026-10-17/tc_energy.json`
  and `comparison.json` as the fixed-code re-run).  Never its `lemma_d` block,
  which is the midpoint record.
- `s4-energy.md` §4–§5 still say "NOT yet certified" and "conditional"; add a
  status banner.

## 8. Cover and provenance (after `a0a1b5ff`, `ccc74f3a`)

**The exact cover (`exact_x_cover`, `exact_gaps`; `s43_parts.py:1509–1715`).**
- *Tube pieces.*
  - The content hash must equal the plan's.
  - The piece must be the step it claims: the drivers' selection, then the
    index, then exact equality of the float labels.  The labels decode at 384
    bits, as `Tube.load` does.
  - The step's radius must be 0.
  - Its interval is [x_k − h_k, x_k] in Fractions.
- *Every cache in use:* all radii 0, and x_{k+1} = x_k − h_k exactly.
- *Sonic pieces:* [−hi, −lo].  *Centre pieces:* [log lo, log hi], with exact
  ordering by ball logs at rising precision.
- *Gaps:* closed-interval sweep with no tolerance.
- *Errors* fail the merge.
- The logic is correct.

**Reproduced.**
- The merge at `75869b91` in my scratch copy gives:
  - 546 substituted, 30 missing, 0 failed, 0 cover errors;
  - gap (−∞, −4.56);
  - per-piece `detail` **identical** to the committed
    `merge_substitute_layer_t14.json` (546 entries, 0 differing fields);
  - `missing`, `failed`, checks and conflicts unchanged.
- `test_s43_merge.py`: 116 passed, 1 skipped.

**Independent check (`cover_independent.py`, no `s43_parts` code).**
- Each piece is matched to its cache step by searching all steps; the match is
  unique for all 546.
- The four caches (344 / 107 / 129 / 60 steps) have radius 0 and exact joins.
- The union is exactly [−5134103575202365/2⁵⁰, 0] = [−4.56, 0], with no gap.
- The smallest cross-region overlaps are
  - 1.80e−16: full-244 / deep-003 (the note says 2.7e−16);
  - 3.12e−5: full-002 / gap-105;
  - 3.79e−5: gap-006 / s50.
- Every `layer_t14` part's cover matches its piece:
  - la, lb;
  - the selection index named in `region_call`;
  - τ ∈ [14, ∞);
  - ε-cells tiling [0, 6].
- `region_call` prints the drivers' **selection** index (`tube_band.steps[0]`
  is cache step 42).  This is cosmetic.

**Provenance.**
- `git_head()` checks return codes and the checkout top level.  Strict
  substitution needs a recorded head and `git_dirty_src is False`.
- `late_modules` is recorded, and fails the write check.
- A static read of the layer path finds no late import: every lazy import is in
  `EAGER`, and `tc_layer_mv` is imported before the snapshot.
- End to end, in the archive copy, I recomputed **sonic/s29** (τ ≥ 14, piece
  norm) and **deep/deep-128** (τ ≥ 14, combine, mv3, 400k/1M):
  - s29: bound −0.4516887835196042, 2515 boxes, 18440 evals;
  - deep-128: bound −0.02165433715776094, 25009 boxes, 199148 evals;
  - both bit-identical to the committed parts, ladder attempts included;
  - both new records have `git_head: null`, `git_dirty_src: null` and
    `late_modules: []`.
- The merge then **rejects** both ("git_head not recorded (not written from a
  checkout)", "git_dirty_src not false").  It substitutes 544, and the exact
  cover reports the two exact holes left behind.  This is the intended
  behaviour.
- All 546 committed parts pass the new strict rule: each records a git head
  and a clean tree.

**Overlap with the centre.**
- c29's hi = YCUT gives x = −4.5 − 3.7e−17.  The deep end is −4.56, so the
  overlap is 0.06.  The merge already computes centre ends as exact logs.
- *Pre-flight for the centre (§4):*
  - the lowest centre layers of dT and dB vanish only at κ ∈ {0, 1};
  - U6's t¹⁰ layer has no κ-root in the strip.
- So nothing structural stops a closing piece at y = 0 from certifying
  dT, dB, nB ≠ 0.

**Cosmetic leftovers (8e).**
- `scope.xi_open_max = Infinity` and `scope.full_strip_tau = [14, 224]` still
  mislabel the open range and the τ range of the substituted pieces.
  - *Failure scenario:* a reader takes "ξ open up to ∞" or "τ only to 224" at
    face value and misstates P4-9's scope.
  - *Fix:* drop `xi_open_max` when the gap is the centre, and report the τ range
    per carrier ("substituted: [14, ∞)").
- `run_centre(..., ncells=6)` still passes an ignored argument.

## 9. Computations

All in the scratch copy of `75869b91`, at nice 10.

| Script (`problems/P4/results/theorem_c/review_s45/`) | What | Result | Time |
|---|---|---|---|
| `check_s45_identities.py` | sympy I1–I7 | all residuals 0; both controls nonzero; I3c numerically spot-checked | < 1 s |
| `g1_order0.py` → `g1_order0.json` | order-0 sonic solve, σ(κ), rank of E | κ\* = −0.0990290820550 ± 2.9e−14 = σ-root; E rank one exactly; det D = −9.6021 | 10 s |
| `centre_pattern.py`, `g1_centre.py` (`.txt`) | centre t = 0 structure | exact zero pattern; det(nP₀ − G₀(κ)) κ-free, identical at five κ | 15 s |
| `centre_layers.py` → `.json` | exact centre layers of TN1r, denT, U6, denU, tt | T + 3 = O(t²), B0 = O(t²), tt = O(t²); denT, denU lowest layers ∝ κ(κ−1)⁴, κ(κ−1)⁶ | 6 s |
| `centre_U6_roots.py` → `.json` | κ-roots of U6's lowest centre layer (float) | none in the strip or Re κ ≥ 6 | 10 s |
| `sonic_residues.py` → `.json` | R_T, R_B, Res tt at x = 0 (balls) | identities hold coefficientwise; γ_min = 1.061671782 | 30 s |
| `factor_checks.py` → `.json` | exact factor structure | all true | 7 s |
| `s42_join.py` → `.txt` | S4-2 join at −4.5 | overlap 4.6e−17, no gap | 2 s |
| `cover_independent.py` → `.txt` | strip x-cover, independent | [−4.56, 0], no gap; parts' covers match | 3 s |
| `s43_parts.py merge --substitute-layer-t14` | the exact merge | as §8; per-piece identical | 1 s |
| `pytest test_s43_merge.py`; the two fixed `test_tc_commutator` tests | regression | 116 passed / 1 skipped; 2 passed | 70 s |
| `s44_layer.py work --force` on s29 and deep-128 | recompute + provenance | bit-identical; rejected by the strict rule outside a checkout | 330 s |

Also checked: the sha256 of every input the assembly lists (I1–I7: results,
cache records, the two part-directory aggregates `e9835a41…` and `974f44e1…`,
and the listed source files) equals the file at HEAD.  The cache content
hashes are re-derived by the merge's cover check.
`git log fda58897..HEAD -- problems/P4/src` is empty.

## 10. Wording fixes (no change of substance)

1. Assembly S8:
   - replace the closing-box and u-hat citations by the residue identities and
     the rigidity lemma of §4;
   - replace "backed exactly by E2" by the centre-layer identities.
2. Assembly A5 / G1: add the centre half (certificate (C)) and the §2 lemma.
3. Class wording in P4-8, P4-9 and P4-10: state the class explicitly (§11).
4. S13 and `s2-theorem-b.md:535`: use κ\* for the complex conjugate.
5. `s4-energy.md`:
   - §2: the xT remainder stored is 1.75e−5, valid on Re κ ≥ 6;
   - §3: "QN/(P³·…)" should read Dp³;
   - §4–§5: add a status banner (Lemma D certified, see
     `lemma_d_certified.json`).
6. Assembly §2: −F_N ranges over [0.3549, 1), not [κ̄, 1); its dip goes below κ̄.
7. Assembly S10: the full-244 / deep-003 overlap is 1.80e−16.
8. Assembly I5: the replay file has pre-ball-fix hashes; cite it as historical.
9. Merge output: see §8, cosmetic leftovers.

## 11. Ledger rows (text only; CLAIMS.md not edited)

"Class" in all three rows:
- nonzero solutions p of the 4D linearised KHA system about the certified EC
  background;
- sonic-point gauge N_p(0) = 0, time dependence e^{κs};
- analytic at the sonic point;
- on the linearised momentum-constraint surface;
- regular at the centre (analytic in t = eˣ with bounded scaled variables).

For Re κ > −0.0990 this equals P4-6's normalized class.

**P4-8 (update).**
- *Claim.* Direct-energy exclusion (Theorem C, part 1): in the class of P4-6
  there is no eigenvalue with Re κ ≥ 6.
- *Method.*
  - The GHJS-transfer Liouville normal form v″ = 𝒬v, with exactly κ-rational 𝒬
    (κ²-density (μ₊ − μ₋)²/4 > 0).
  - Boundary terms vanish in the class.  At the sonic point, xT → 1 + σ(κ) with
    Re σ < 0.  At the centre, T = −3 + O(t²) exactly.
  - Conditions (a)–(d) are certified pointwise on all of x ∈ (−∞, 0):
    - order-2 jet arithmetic on [−4.5 − 8.3×10⁻¹⁷, 0]: tube 14520 + sonic
      ladder 16666 + the closing box, the sonic point included; exact
      identities QNm = D¹⁰·QN_red and QN_red(u₀) ≡ 0 in ℚ(V₀, √3);
    - Lemma D on y = eˣ ∈ (0, YCUT], YCUT = 0.011108996538242306 < e^{−4.5}:
      Taylor models in y on the certified (a\*, μ\*) balls, 1100/1100
      subranges; exact layer identity E2;
    - the two regions overlap by 4.6×10⁻¹⁷.
- *E-clause.*
  - On [6, 15] × [−14, 14], E(κ) ≠ 0 follows by P4-6.
  - Replace by "hence E(κ) ≠ 0 on {Re κ ≥ 6}" once the G1 lemma (review §2) is
    committed.
- *Level.* **CERTIFIED**, with the Lemma D condition dropped, once this
  review's checks `centre_layers` and `s42_join` are committed.  Until then the
  row keeps its current status.
- *Evidence.*
  - `problems/P4/notes/s4-energy.md` §1–§3, §6;
  - `problems/P4/results/theorem_c/parts/lemma_d_certified.json`;
  - `.../parts/checks.json`;
  - `.../tc_energy.json`: the tube, sonic and closing blocks only; its
    `lemma_d` block is the superseded midpoint record;
  - `problems/P4/results/recert_2026-10-17/{tc_energy,comparison}.json`;
  - this review (suggested path `docs/review-2026-10-17-p4-s4-5-adversarial.md`)
    and its committed checks.

**P4-9 (new).**
- *Claim.* Strip exclusion (Theorem C, part 2): in the class of P4-8 there is no
  eigenvalue with 0 ≤ Re κ ≤ 6 and |Im κ| ≥ 14.
- *Method.*
  - A once-commuted (m = 1) energy identity for χ₁ = χ′, with weight
    g = exp∫(Re A1 − c·tt) and c = −1/32.
  - The pointwise condition S = K + ν²/4 < 0, with dT, dB, nB ≠ 0, is
    certified on piece × ε ∈ [0, 6] × τ ∈ [14, ∞) for every piece of an exact
    cover of x ∈ (−∞, 0).
  - Boundary terms: γ = 1 + σ₀ + (1 − σ₀)(ε + c) ≥ 1.0617 at the sonic point;
    BT = O(t^{3+k}) → 0 at the centre, with B0 = O(t^k), k ≥ 2.
  - Im κ ≤ −14 by conjugation.
- *Level now.* **CONJECTURED.**
  - 546 of 576 pieces are certified: x ∈ [−4.56, 0], exact cover, S4-4
    `layer_t14` parts under `merge --substitute-layer-t14`.
  - The 30 centre pieces (x ≤ −4.5) are missing.
- *Target.* CERTIFIED when the centre pieces pass on the current code with the
  certified centre balls and the merge reports `strip_exclusion_certified =
  true`.  This review accepts S5–S8, S12 and S13 as routine, given its committed
  checks.
- *Evidence.*
  - notes: `s4-commutator.md` §2–§3, `s4-4-layer-engine.md` §1–§5 and §12–§14,
    `s4-5-assembly.md` §4.1;
  - `parts/layer_t14/` (aggregate `e9835a41…`);
  - `layer_recompute/merge_substitute_layer_t14.json`;
  - `parts/checks.json`;
  - the two 2026-10-17 review documents and this review's checks.

**P4-10 (new).**
- *Claim.* Theorem C (eigenvalue form) and the mode count.
  - In the class of P4-8 there is no eigenvalue in
    Ω = {Re κ ≥ 0} \ [0, 15] × [−14, 14].  P4-8 covers Re κ ≥ 6; P4-9 and
    conjugation cover the rest.
  - With P4-6, the eigenvalues with Re κ ≥ 0 are exactly
    κ̄ ∈ [0.35569919907, 0.35569920835] (pure gauge) and
    κ₁ ∈ [2.81055254439, 2.81055255326].  Both are simple zeros of E.
  - So κ₁ is the unique non-gauge mode with Re κ ≥ 0, and
    γ = 1/κ₁ ∈ [0.3558019218, 0.3558019231].
  - Not claimed: kink modes, non-spherical modes, nonlinear stability,
    Re κ < 0.
  - The E-form "E(κ) ≠ 0 on Ω" is a corollary once the G1 lemma is committed.
- *Level now.* **CONJECTURED** (it inherits P4-9).
- *Target.* CERTIFIED with P4-9.  If P4-8 or P4-9 ends at WRITTEN_PROOF, so does
  P4-10.
- *Evidence.* The P4-6, P4-8 and P4-9 evidence; `s4-5-assembly.md`; this
  review.

**Policy point (unchanged from assembly §5 item 15).**  README defines CERTIFIED
by a certificate "committed under `certificates/`".  P4 keeps its artifacts
under `results/`, as the existing CERTIFIED rows P4-3 … P4-6 already do.  The
owner should note or change this.
