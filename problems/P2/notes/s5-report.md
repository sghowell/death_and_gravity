# P2 S5 report: the assembly P2-4 and the conversion P2-5 at P_B

Status: S5 exit, 2026-10-26. Branch `p2/s5-assembly`, from main at `3dd0b2bc`.
- Plan: `notes/s5-plan.md`, committed first at `e2c75aa6`, with no amendments.
- Written proof and citation analysis: `notes/s5-assembly.md`.
- This report has not yet been independently reviewed. S5's adversarial review (FORMULATION §6) is
  still to come.

## 1. Verdict

- **P2-4 at P_B holds, at WRITTEN_PROOF.**
  - For every l, W_l has no zero in S_δ = {Im ω ≥ −0.55κ₋} other than the simple zero of W₀ at 0.
  - Equivalently, the forward problem on Ω₂₃ has no resonance in S_δ∖{0}, and 0 is a resonance of
    multiplicity 1 with constant resonant states.
- **β at P_B: 0.55 ≤ β ≤ 0.90938098.** The bound is two-sided but not sharp. P2-6 stays the stretch
  goal.
- **P2-5 at P_B closes by citation, at WRITTEN_PROOF, conversion cited.** It is not conditional:
  - [HV17]'s α hypothesis is about the forward problem on Ω₂₃. This is quoted from §2.5 of the TeX,
    and the next sentence there excludes the global P̂(σ)⁻¹ in so many words.
  - P2-D1′ identifies that notion with the zeros of W_l.
  - The global resonances are removed in r > r₁ by [HV17] Lemma 2.14(1), inside [HV17]'s own proof.
- The plan's stop rule is not triggered.

## 2. Commits (branch `p2/s5-assembly`)

| commit | content |
|---|---|
| `e2c75aa6` | `notes/s5-plan.md`: deliverables, the assembly-check design and its negative controls, the hypotheses to match, acceptance criteria, the stop rule |
| `7aa9fbd0` | `src/p2/assembly.py` (`python -m p2.assembly [verify]`), `tests/test_assembly.py` |
| `942f89c1` | `certificates/assembly/assembly_PB.json`, computed at `7aa9fbd0` on a clean tree |
| `a04a4323` | `notes/s5-assembly.md` (the P2-4 proof and the P2-5 analysis with quotes); FORMULATION v1.0 |
| this commit | `notes/s5-report.md` |

Only `problems/P2/` was changed. Not touched: `CLAIMS.md`, the .tex, main.

## 3. The assembly check

`python -m p2.assembly` reads 54 committed certificates and pins each by its sha256:
- params;
- lg and winding for l = 0–19;
- the part-1 cover and the regime certificates;
- the P2-0 identities, Z4, and the Krawczyk/KAT-6 enclosures.

It uses exact rationals and 256-bit Arb. **All 366 checks pass**, in 2.7 s of CPU. The checks:
- **l ≤ 19.**
  - The a_l of R′_l is identical to P2-2's stored threshold.
  - b − d_δ = 4.546·10⁻¹⁰ > 0, e > 0, and b < κ₊.
  - The 2351 stored discs trace the right half of ∂R′_l exactly and clear E₊ ∪ E_c.
  - The re-summed windings are 1 at l = 0 and 0 otherwise.
  - The lg criteria are below 2.
- **l ≥ 20.**
  - `part1.cover()` replays identically.
  - Every regime certificate passes, with h₀ = 2/41.
  - l₀ = ⌈1/h₀ − ½⌉ = 20, so no l is missed.
  - L ∪ S ∪ T ∪ A covers Re ω ∈ [0, ∞) for every h ≤ h₀, with S's top at h₀ at 90.75κ₋ against
    L's 8κ₋.
  - The Im range is [−d_δ, 0].
- **Poles.** d_δ < κ₊ (margin 0.01161) and d_δ < κ_c (margin 0.06998).
- **The zero mode.** 0 ∈ int R′₀, and W₀′(0) = 5.5491372i.
- **β.** β ≤ 0.90938098, with margin 1.005·10⁻⁹.
- **The HV17 and CF17 numerics** (§4).

**The negative controls N1–N8 all fail as required.** The tightest is N1: b moved to 6280878/2²⁸,
just below d_δ, fails on `C2.l0.b_ge_d` and nothing else. The others:
- N2: `winding_PB_l7.json` removed;
- N3: a₅ − 2⁻²⁴;
- N4: h₀ = 2/43, so l₀ = 21;
- N5: e = 0;
- N6: a wrong winding;
- N7: a failed regime;
- N8: γ₀ and κ₊ lowered, so that HV17's α < γ₀ and CF17's κ₋ < 2κ₊ are exercised.

**Replays in S5.**
- `python -m p2.assembly verify`: identical.
- `python -m p2.params verify`: OK.
- `python -m p2.lg_bounds verify PB 19`: OK.
- The P2 fast suite: 90 passed, 17 skipped. That is 86 + 4 new tests.

## 4. [HV17] and [CF17]: the match

Full quotes, with labels and TeX line numbers, are in `s5-assembly.md` §2.

| hypothesis | source | at P_B |
|---|---|---|
| H1: non-degenerate RNdS, Q ≠ 0, dim 4 | Thm 1.1 [ThmIntroMain]; Def 2.1 [DefRNdSNonDegenerate] | μ = F; Λ, M, Q > 0; r_n < 0 < r₋ < r₊ < r_c, all simple (P2-par) |
| H2: □_g u = 0, smooth data on H_I | Thm 1.1; (EqRNdSMfdHyper) | massless scalar ✓. Data on the one-ended H_I; conclusions in its domain of dependence. Two-ended data: near every point of CH⁺ except the bifurcation sphere (Killing translates); up to it, [HV17] only sketches it, and it is not claimed |
| **H3: resonance notion of "0 is the only resonance in Im σ ≥ −α"** | §2.5 [SubsecRNdSAsymp]; Lemma 2.14 [LemmaRNdSResonanceSupp]; proof of Prop 2.17 | **Forward problem on Ω₂₃, the localized notion = P2-D1′.** The same paragraph: "(This does *not* mean that the global problem for P does not have other resonances in this half space!)". Global resonances that are not Ω₂₃ resonances are removed in r > r₁ by Lemma 2.14(1) |
| H4: 0 simple, constant resonant states | §2.5; Lemma 2.14(2) | [HV17]'s own claim, which rests on references for Schwarzschild–de Sitter and Kerr–de Sitter. Independently: D4 (m₀(0) = 1), P2-0 and P2-D1′ |
| H5: α < γ₀ | §2.5; Prop 2.12 [PropRNdSSemiclassical2]; §2.7 | γ₀ is D6's. d_δ/γ₀ = 0.327 (margin 0.0482, P2-par) |
| H6: exponent | (EqRNdSRadialPointBeta), (EqRNdSSurfGrav); Prop 2.17; Thm 2.24 | β₁ = 1/κ₋, so ½ + α/κ₋ = 21/20 exactly: H^{1.05−0} ⊂ H¹ |
| H7: implicit points | proof of Prop 2.17 | The contour Im σ = −α may meet global poles. They are harmless (Lemma 2.14(1), or run at α̂ ↑ α). This is a remark on [HV17]'s proof, not a missing hypothesis |
| F1: subextremal RNdS, φ ∈ C^∞(M_int ∪ H⁺) | [CF17] §2.1; Thm 3.1 [mainThm] | D = F, e = Q ✓. Data on a Cauchy surface Σ |
| F2: Price law (Assumption 2.2 [pricesAss]), rate p | — | **p = α**, from Thm 1.1's bounds on stationary derivatives tangent to CH⁺, using T = ∂_v on H⁺, t₀ − v bounded near H⁺, and abs(∇̸ψ)² = r⁻²Σ(Ω_iψ)² |
| F3: 2 min{p, κ₊} > κ₋ | Thm 3.1 (mainCond) | ⟺ 2p > κ₋ and 2κ₊ > κ₋. 2p/κ₋ = 11/10 is the **tautology** (any δ > 0). 2κ₊/κ₋ = 1.646 is the substantive part (P2-par). The CF17 route's content is therefore the Price law with p > κ₋/2, i.e. P2-4 |
| F4: κ conventions | (surfGrav) | κ_* = ½ abs(D′(r_*)), the same as ours |

## 5. Levels

| row | level | why |
|---|---|---|
| P2-4 | **WRITTEN_PROOF** | The minimum of P2-par (CERTIFIED), P2-0, P2-D1′, P2-1 (WRITTEN_PROOF), P2-2 and P2-3B (CERTIFIED). See below |
| P2-5 | **WRITTEN_PROOF, conversion cited** | P2-4's level. [HV17] and [CF17] are cited, not re-proven |
| β ≥ 0.55 | WRITTEN_PROOF | P2-4 |
| β ≤ 0.90938098 | P2-D1′'s level (WRITTEN_PROOF) | P2-K6, unchanged |

Why P2-4's level is the minimum:
- Each point of S_δ∖{0} must be excluded for every l.
- The pieces doing that are all WRITTEN_PROOF:
  - on U, P2-0 for every l;
  - P2-1 for the infinitely many l ≥ 20;
  - the resonance form, P2-D1′.
- No CERTIFIED row replaces any of them.
- The assembly's own steps are routine or machine-checked.

For l ≤ 19 the lower strip is CERTIFIED. P2-4 would be CERTIFIED only if P2-0's written steps,
P2-1's W1–W2 and D1′ were machine-checked.

## 6. Ledger text: paste-ready, not applied

Fill in the ⟨…⟩ fields after S5's independent review. The column order is
`| id | problem | statement | level | evidence | updated |`.

```
| P2-4 | P2 | At P_B (ΛM² = 3/50, Q = (999/1000)Q_max, δ = 1/20), for every l ≥ 0: W_l has no zero in S_δ = {Im ω ≥ −(½+δ)κ₋} other than the zero of W₀ at ω = 0, which is simple (m₀(0) = 1). Equivalently (P2-D1′), the forward problem on Ω₂₃ (HV17 §2.5) has no resonance in S_δ∖{0}, and ω = 0 is a resonance of multiplicity 1 whose resonant states are the constants. Hence D3's gap satisfies β ≥ ½ + δ = 0.55; with P2-K6, 0.55 ≤ β ≤ 0.90938098 (two-sided, not sharp; the two-sided enclosure is the stretch goal P2-6). Assembled from P2-0 (Im ω > 0, every l), the conjugation symmetry (Re ω < 0), P2-1 (l ≥ 20), P2-2 with P2-3B (l ≤ 19: abs(Re ω) ≥ a_l, and R′_l) and P2-par (E₊ ∪ E_c below S_δ; S_δ inside P2-D1′'s multiplicity window). Not P_A | WRITTEN_PROOF, the minimum of the ingredients: P2-0, P2-D1′ and P2-1 are WRITTEN_PROOF; P2-par, P2-2 and P2-3B are CERTIFIED (for l ≤ 19 the lower strip −(½+δ)κ₋ ≤ Im ω ≤ 0 is CERTIFIED). Machine-checked bookkeeping (exact rationals and Arb, 256-bit; 366 checks; the 54 input certificates pinned by sha256): for each l ≤ 19, the a_l of R′_l is identical to P2-2's stored threshold, b − (½+δ)κ₋ = 4.546e−10 > 0, e > 0 and b < κ₊; the 2351 stored discs trace the right half of ∂R′_l exactly and clear E₊ ∪ E_c; the re-summed windings are 1 (l = 0) and 0. P2-1's cover replays identically; h₀ = 2/41 gives l₀ = 20, so no l is missed; the four regimes cover Re ω ≥ 0 for every h ≤ h₀. (½+δ)κ₋ < κ₊, κ_c. Negative controls N1–N8 fail as required: b just below (½+δ)κ₋; one l's certificate removed; a₅ lowered by 2⁻²⁴; h₀ = 2/43; e = 0; a wrong winding; a failed regime; γ₀ and κ₊ lowered. Written (notes/s5-assembly.md §1, routine): the decomposition of S_δ, the argument principle given the certified counts, the D3 infimum. Independent review ⟨date⟩: ⟨verdict⟩ | problems/P2/src/p2/assembly.py; problems/P2/certificates/assembly/assembly_PB.json; problems/P2/notes/s5-plan.md, s5-assembly.md, s5-report.md; problems/P2/tests/test_assembly.py; ⟨review document⟩ | ⟨date⟩ |
| P2-5 | P2 | At P_B, with α = (½+δ)κ₋ = 0.55κ₋ (< γ₀ = 1.683κ₋), by citation. (i) [HV17] (arXiv:1512.08004) Thm 1.1 in its quantitative form (Prop 2.17, Cor 2.19, Thm 2.24): every solution of □_g u = 0 with smooth initial data on HV17's Cauchy surface H_I equals u₀ + u′ in the domain of dependence of H_I, with u₀ ∈ ℂ and abs(u′) ≤ Ce^{−αt₀} uniformly in r > r₋ (also after any number of stationary derivatives tangent to the Cauchy horizon). u is continuous up to the Cauchy horizon, and u′ ∈ e^{−αt₀}H^{1/2+α/κ₋−0} = e^{−αt₀}H^{1.05−0} in t₀ > 0; in particular u ∈ H¹_loc up to the Cauchy horizon. (ii) [CF17] (arXiv:1607.01018) Thm 3.1 with p = α, for smooth data on a Cauchy surface: by (i) the Price law (Assumption 2.2) holds with p = α, and 2 min{p, κ₊} = 1.1κ₋ > κ₋; so the energy is bounded up to the Cauchy horizon, and the solution extends to it in C⁰ ∩ H¹_loc. This is the linear analogue at P_B of the failure of Christodoulou's formulation of SCC for smooth data. Massless scalar on fixed RNdS only; for two-ended data, not uniformly up to the bifurcation sphere; not P_A | WRITTEN_PROOF (P2-4's level), conversion cited: neither theorem is re-proven. HV17's α hypothesis (§2.5), quoted from the TeX: "The only resonance of the forward problem in Ω₂₃ in Im σ ≥ 0 is a simple resonance at σ = 0 … and there exists α > 0 such that 0 is the only resonance in Im σ ≥ −α. (This does not mean that the global problem for P does not have other resonances in this half space!)" The global resonances are removed in r > r₁ by HV17 Lemma 2.14(1), in the proof of Prop 2.17. So P2-4 with P2-D1′ supplies the hypothesis exactly, with no extra argument. 0 simple with constant states is HV17 Lemma 2.14(2), and independently D4 with P2-0 and P2-D1′. From P2-par: α < γ₀ (margin 0.0482) and κ₋ < 2κ₊ (margin 0.0275). The exponent ½ + α/κ₋ = 21/20 is exact. In CF17's condition, 2p > κ₋ is the tautology noted in P2-par's review, so the CF17 route's content is the Price law with p > κ₋/2, i.e. P2-4. The numerical hypotheses are checked in certificates/assembly/assembly_PB.json (C7, C8, with negative controls). Independent review ⟨date⟩: ⟨verdict⟩ | problems/P2/notes/s5-assembly.md §2; problems/P2/certificates/assembly/assembly_PB.json; ⟨review document⟩ | ⟨date⟩ |
```

Optional, for the user: the P2-1 row's parenthesis "(the assembly P2-4 is S5)" could later read
"(assembled in P2-4)".

## 7. Proposed P2 status entry for the .tex: text only, not applied

This is a revision-log item in the style of v1.2, a new bibliography entry, and an optional sentence
for the P2 Status paragraph. Fill in the ⟨…⟩ fields after the review and the approval.

```latex
\item[v1.3] ⟨Month day⟩, 2026 (proposed October 26, 2026 by the P2 stage-S5 assembly, \texttt{problems/\allowbreak P2/\allowbreak notes/\allowbreak s5-report.md}; approved ⟨date⟩). P2 only; the other problems are unchanged. \emph{P2 Status.} (a) At the primary pre-registered point $\Lambda M^2=3/50$, $Q/Q_{\max}=999/1000$, the resonance-free region of the entry point holds with $\delta=1/20$ for the linear massless scalar field: for every $\ell$, the strip $\{\operatorname{Im}\omega\ge-0.55\,\kappa_-\}\setminus\{0\}$ contains no quasinormal mode, and $\omega=0$ is a simple resonance with constant resonant states ($0.55\,\kappa_-\approx0.33\,\gamma_0$). Hence $0.55\le\beta\le0.90938098$ there; the upper bound comes from an interval-certified $\ell=0$ near-extremal mode. (b) The proof combines a written mode-stability argument on the closed upper half-plane; a large-$\ell$ exclusion, uniform in $\ell\ge20$, whose inequalities are certified in interval arithmetic and whose analytic steps are written; certified large-$|\operatorname{Re}\omega|$ exclusions and validated argument-principle counts for $\ell\le19$; and a written identification of the zeros of the radial Wronskian with the resonances used in \cite{hintzvasy2017}. The assembly is machine-checked. The result is a written proof with machine-checked certificates, independently reviewed ⟨date⟩; it is not a fully machine-checked certificate. (c) By \cite[Thm.~1.1 and \S\S2.5--2.6]{hintzvasy2017}, whose spectral hypothesis concerns the forward problem near the black-hole exterior (the notion established in (a)), every solution of the linear wave equation with smooth Cauchy data on the one-ended Cauchy surface of \cite{hintzvasy2017} lies in $H^{1.05-0}$, in particular in $H^1_{\mathrm{loc}}$, at the Cauchy horizon. \cite[Thm.~3.1]{costafranzen2017} gives bounded energy and a $C^0\cap H^1_{\mathrm{loc}}$ extension. This is the linear analogue of the failure of the Christodoulou formulation for smooth data at this point. (d) Not established: the Einstein--Maxwell problem; the second point $Q/Q_{\max}=249/250$; a two-sided enclosure of $\beta$ (the stretch goal); any statement at other parameters, or uniform up to the bifurcation sphere of the Cauchy horizon. \emph{Bibliography.} Added \cite{costafranzen2017}.

\bibitem{costafranzen2017} J.~L.~Costa and A.~T.~Franzen, ``Bounded energy waves on the black hole interior of Reissner--Nordstr\"om--de Sitter,'' Ann. Henri Poincar\'e \textbf{18}, 3371 (2017), arXiv:1607.01018.
```

Optional sentence for the end of the P2 Status paragraph, after "…yields only
$\beta\ge1/2$.":

```latex
At the primary entry-point parameters this region has since been established with $\delta=1/20$ (revision v1.3).
```

## 8. What S5 does not establish

- **Nothing at P_A** (P2-7). Part (i) has not been run there.
- **No two-sided enclosure of β to 10⁻⁸** (P2-6, the stretch goal). The upper bound is one certified
  mode. There is no strict β > 0.55, and no claim about the spectrum between 0.55κ₋ and 0.909κ₋
  beyond P2-4.
- **No CERTIFIED level** for P2-4 or P2-5. P2-0's written steps, P2-1's W1–W2 and D1′ remain prose.
- **No re-proof of [HV17] or [CF17].** In particular, [HV17] Lemma 2.14(1)'s general clause and the
  order-1 claim at 0 are cited; the latter is also P2-0/D4 with D1′ at P_B.
- **For two-ended data, nothing uniform up to the bifurcation sphere** of the Cauchy horizon.
  [HV17]'s "minor extension" is not used.
- **The scope fences of FORMULATION §7 are unchanged:**
  - no nonlinear Einstein–Maxwell statement, and the [CGNS18]/[Ros25] hypotheses are not
    established;
  - no other fields or spacetimes;
  - no sharpness: nothing about generic solutions being no smoother than H^{1/2+β};
  - no other parameters or neighbourhoods;
  - nothing on [Hin25] Conj. 1.5.
- **No independent review yet.** The ledger rows and the .tex text are proposals.

## 9. For the reviewer: what to attack first

1. **The resonance notion** (`s5-assembly.md` §2.3). Attack:
   - the reading of [HV17] §2.5;
   - whether the proof of Prop 2.17 needs the global notion anywhere else;
   - the reliance on Lemma 2.14(1)'s "more general" clause, which [HV17] proves only in outline.

   This is where P2-5 would break.
2. **The P2-4 seams**:
   - the identity of the a_l thresholds;
   - b against d_δ, with an absolute margin of 4.5·10⁻¹⁰;
   - the line Im ω = 0, held by both P2-0 and the lower-strip certificates;
   - the boundary between l = 19 and l = 20, at h₀ = 2/41 exactly;
   - 0 ∈ int R′₀;
   - whether the four regime ranges really cover the closed half-strip for every h ≤ h₀.
3. **The Price-law derivation (F2)**:
   - t₀ against v near H⁺;
   - T = ∂_v on H⁺;
   - which component of H⁺ [HV17]'s H_I crosses;
   - the compact initial piece.
4. **The setting (H2).** Data on H_I against data on a complete Cauchy surface, and the bifurcation
   sphere.
5. **The contour remark (H7)**, and the claim that the level is the minimum of the ingredients.

## 10. Budget

The assembly runs about 3 s of CPU each (about 6 runs). Each fast-suite run takes 37 s, and the
replays about 4 s. The total is under 0.05 core-h, at nice 10, with one worker.
