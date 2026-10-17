# P2 S0 gate: go/no-go memo for the certified resonance-free strip

Status: S0 memo, 2026-10-16. Nothing here is certified. Numbers without a ball
radius are floats from `src/s0_sanity.py` (mpmath, 25–40 digits, about 5 CPU-min
in total, one core). Inputs: `../FORMULATION.md` v0,
`literature-digest.md` (bracketed keys), P4 Theorem B and C notes
(`problems/P4/notes/s2-theorem-b.md`, `s4-energy.md`), and the commutator note
(`git show p4/theorem-c-s4-3:problems/P4/notes/s4-commutator.md`).

## 0. Verdict

**GO, with a modified entry point** (proposed v1.2 in §7; needs your approval
before S1 freezes anything that depends on it).

Part (ii) (large |Re ω| at fixed l) and part (iii) (the validated count) use
technology this repo already has. Part (iii) is a simpler ODE than P4's: the
coefficients are exact algebraic numbers and no numerical background tube is
needed.

Part (i) (uniform large l) is the research core. No explicit-constant
barrier-top result exists for a *complex* spectral parameter. It is still a
one-dimensional problem with an explicit analytic potential, and how hard it is
depends almost entirely on the parameter point:
- At the point in [CCDHJ18] Table I (P_A), β ≈ 0.764 is attained only in the
  eikonal limit l → ∞. Part (i) must then exclude a strip of width 0.67 of the
  trapping gap γ₀.
- In the near-extremal window of [CCDHJ18] Fig. 3 (P_B, Q/Q_max = 0.999),
  β ≈ 0.909 is attained by one isolated l = 0 mode. Part (i) then only needs
  0.33·γ₀, a threefold margin.

The modification pre-registers P_B as the primary point. It also corrects the
strip statement (δ > 0) and requires part (i) to hold over every regime of
Re ω, not just the barrier top.

A kill gate at S2 (§7) stops the program before heavy investment if the
barrier-top certificate cannot be projected to l₀ ≤ 100.

## 1. What must be proven

At a frozen point X (FORMULATION §1.3): no zero of W_l in
S_δ = {Im ω ≥ −(½+δ)κ₋}∖{0} for any l, and ω = 0 a simple zero of W₀. Then
[HV17] gives H^{1+δ−0} at the Cauchy horizon (P2-5). The closed upper
half-plane is classical (§2). The work is the strip −(½+δ)κ₋ ≤ Im ω < 0,
split into (i) l ≥ l₀, (ii) l < l₀ with |Re ω| ≥ W_l, and (iii) the compact
boxes R_l.

Decisive numbers (float):

| | P_B (999/1000), primary | P_A (249/250), second |
|---|---|---|
| dominant mode, β | l=0 NE, −0.909381i κ₋: **β ≈ 0.9094** | PS as l→∞: **β = γ₀/κ₋ = 0.76430** |
| next shallowest | l=0 NE n=1 (1.855); l=1 NE (1.912); PS (≥ 1.683) | Table I l=100: 0.76431; l=10: 0.76492 |
| strip edge, δ | 0.55 κ₋ (δ = 1/20) | 0.51 κ₋ (δ = 1/100) |
| margin to nearest QNM | 0.36 κ₋ | 0.25 κ₋ (only in the limit l → ∞) |
| **part (i) width / γ₀** | **0.327** (0.547 for two-sided β, P2-6) | **0.667** |
| κ₊/κ₋, κ_c/κ₋ | 0.823, 2.195 | 0.677, 0.990 |
| H^s needed at r₊ / r_c (hyperboloidal) | s > 1.17 / s > 0.75 | s > 1.25 / s > 1.02 |
| Jost margin κ₊ − (½+δ)κ₋ | 0.273 κ₋ | 0.167 κ₋ |

The Frobenius exceptional sets E₊ = {−inκ₊} and E_c = {−inκ_c} lie outside
S_δ at both points.

## 2. Why the damped strip is hard (and what P4 does and does not give)

In hyperboloidal variables (τ = t − h(r), Fh′ = η, η(r₊) = −1, η(r_c) = +1,
R = e^{iωh}u), the radial ODE becomes

    (r²F u′)′ + iω[2r²η u′ + (r²η)′u] + ω² r²(1−η²)F⁻¹ u − L u = 0,

with u analytic at both horizons for a QNM. Multiply by ū and integrate over
(r₊, r_c). The boundary terms of the second-order part vanish because F = 0
there. With
- D = ∫r²F|u′|², N = ∫|u|², 𝒲 = ∫r²(1−η²)F⁻¹|u|² > 0,
- X = ∫r²η Im(u′ū),
- B = r_c²|u(r_c)|² + r₊²|u(r₊)|² ≥ 0,
- ω = a + ib,

the real and imaginary parts are

    −D − bB − 2aX + (a²−b²)𝒲 − LN = 0,     aB − 2bX + 2ab𝒲 = 0,

so that, for b ≠ 0,

    **D + L·N + |ω|²𝒲 = −(|ω|²/b)·B.**                       (E)

- **b > 0.** Both sides of (E) must vanish, so u ≡ 0 (or (l,ω) = (0,0) with u
  constant).
- **b = 0, a ≠ 0.** The imaginary part gives B = 0, so u(r₊) = 0. This
  contradicts the normalization of the analytic Frobenius branch, so u ≡ 0.

This is P2-0: exact, no computation. It is the T-energy flux balance
(dE/dτ = −flux through the horizons). At S0, two pieces were checked as exact
identities by a sympy scratch run: the transformed ODE (for generic F, h, u)
and the elimination of X that gives (E). The integration by parts is the
written step. S1 turns this into the committed P2-0 chain.

**For b < 0, (E) is consistent with any positive left side.** The flux term
(|ω|²/|b|)B pays for the bulk. So energy identities exclude growing modes
(P4 S4-2's situation, Re κ ≥ 0, where the boundary terms vanish) but not
damped ones. Three structural facts make the damped strip harder than P4's
exclusion:
1. **Horizon thresholds.** QNMs are H^s eigenfunctions only for
   Im ω > −(s−½)κ at each horizon ([Vas13], [War15]). H¹ is below threshold
   at r₊ at both points, and also at r_c for P_A (table above). Any multiplier
   method must work at the H² level, which means commuting once with the
   redshift field. This is the analogue of P4's "commute once" step
   (s4-commutator §2), so it is an engineering cost, not a blocker.
2. **Trapping.** At the photon sphere a real Morawetz multiplier proves only an
   l-dependent strip of width O(1/(l log l)) or O(1/l), because the commutator
   is bounded below by about h² at the hyperbolic fixed point. A uniform strip
   of width c·γ₀ needs the inverted-oscillator structure of the trapped set:
   the uncertainty lower bound hλ at the top ([WZ11], [Dya16]). P4's
   large-|κ| problem had no trapping ("single WKB term",
   theorem-c-plan §2.3), so this is genuinely new.
3. **Uniformity in h = 1/(l+½) → 0 with coupled scaling.** At the top the
   natural variable is σ = (r* − r*_P)/h^{1/2}. P4 S4-3's open item (the
   τ > 224 near-sonic sliver, with coupled scaling ξτ where ξ is the distance
   to the sonic point) is the in-repo precedent for exactly this kind of
   corner.

## 3. Part (ii): large |Re ω| at each fixed l < l₀

| method | gives | missing | explicit constants? | size |
|---|---|---|---|---|
| **N1** Born/Volterra (Jost solutions as perturbations of e^{∓iωr*}, matched at r*_P) | W_l(ω) = 2iω(1+ε), with ε ≤ exp(I_l/\|ω\|) − 1 and I_l = ½∫(1+e^{2\|b\|\|x−x_P\|})\|V_l\|dx; converges iff \|b\| < κ₊, κ_c (true, Jost margins above) | the threshold grows like O(l(l+1)), so the part-(iii) boxes grow quadratically in l | yes: integrals of explicit functions of r, Arb quadrature | ≤ 1 session; fallback only |
| **N2** Liouville–Green above the barrier with Olver error control on the real progressive path ([Olv74] Ch. 6 §11) | W_l ≠ 0 for \|Re ω\| ≥ W_l ≈ (1+ε)√max V_l, about (l+½)Ω_P(1+ε); the real line is progressive because Im(ω²−V) = 2ab has fixed sign; the error-control variation is about ∫\|V″\|/\|f\|^{3/2} + ∫\|V′\|²/\|f\|^{5/2} and is finite at both ends for \|b\| < κ₊, κ_c | nothing conceptual; l = 0, 1 need larger ε (weak barrier, negative V₀ near r_c), which is cheap | yes: error-control integrals as Arb quadrature on r-boxes, with endpoint series at r₊ and r_c | about 1 session; **recommended** |
| **N3** hyperboloidal energy and Morawetz identity at high frequency (P4 S4-2 style) | non-trapping exclusion at fixed l | function-space thresholds (H²) and multiplier design | yes | 1–2 sessions; no advantage over N2 |

**Verdict (ii): achievable with explicit constants, low risk, about 1 session.**
Planning value W_l = 1.25(l+½)Ω_P + 2κ₋.

## 4. Part (i): uniform large-l exclusion (the research core)

Semiclassical form: h²Φ″ = (U(x) + h²W̃(x) − z)Φ, with h = 1/(l+½),
U = F/r² (one maximum U₀ at r_P), z = h²ω², and W̃ = FF′/r − U/4. A uniform
proof needs four regimes:
1. **Above the barrier**, Re z > U₀(1+ε): handled by N2 with error O(h).
2. **The barrier top**, a PCF/Weber layer through the coalescence of two
   complex turning points.
3. **Below the barrier**: Airy turning layers, plus forbidden regions where the
   Riccati half-plane {Re y > 0} is forward-invariant whenever
   Re(U − z) > 0 (proved in one line: on Re y = 0, (Re y)′ > 0).
4. **Low frequency**: turning points pushed into the exponential tails. On
   P_B's throat the comparison is with modified Bessel I_{−iω/κ₊}
   (V ≈ c·e^{2κ₊x}).

At the top, in the scaled model, the strip corresponds to
0 ≤ Im μ ≤ c* (c* = the width/γ₀ fraction), while the model resonances sit at
Im μ = 1, 3, 5, …. **P_B leaves a gap of 0.67 (c* = 0.33); P_A leaves 0.33.**

| method | gives | missing | explicit constants? | size, risk |
|---|---|---|---|---|
| **M1** Olver–Dunster uniform parabolic-cylinder approximation ([Olv75]; [Dun21]; [Dun25], arXiv:2504.19405) | uniform approximation through turning-point coalescence with explicit error bounds O(1/u), u = l+½; matching to LG tails | **complex coalescence parameter.** [Olv75] and [Dun25] assume a real parameter (two real or conjugate turning points); for complex ω the turning points are complex and not conjugate, so the error-control theory must be re-derived on suitable progressive paths. The Liouville map ζ(x) is implicit and needs Taylor-model enclosures | yes in principle; the error-control functions are explicit integrals | 4–6 sessions of analysis plus 1–2 of certification; risk **medium** |
| **M2** Riccati invariant tubes (P4-style certified pointwise sign condition in (r, h) plus the spectral parameters): track y₊ = Φ₊′/Φ₊ and y_c in discs around a computed centre ỹ; a disc is invariant if \|ρ̃\| − 2r Re ỹ + r² ≤ r′ pointwise (ρ̃ = Riccati defect of ỹ); QNM ⇔ y₊(x_P) = y_c(x_P) | everything reduces to pointwise inequalities in (x, h, Re ω, Im ω), checkable in Arb boxes with closing boxes at h = 0 in scaled variables; no oscillatory quadrature | a good centre near the top. The pure quadratic model gives relative error O(h^{1/2}σ³) and patches to LG only at O(h^{1/4}) (§4.1); the centre must be M1-quality (uniform, O(h)) or a certified high-order fit. Four parameters; neutral stability (Re ỹ ≈ 0) through the allowed region forces Grönwall growth, finite because \|b\| < κ₊; the projective (Riemann-sphere) form is needed wherever Φ can vanish | yes | 4–8 sessions; risk **medium**; the most "native" route |
| **M3** complex deformation and coercivity: rotate r* by about π/4 at the top and scale the tails (outgoing solutions decay); a variable-phase multiplier reduces exclusion to positivity of a 1D Schrödinger form, certified by a ground-state (Agmon–Allegretto–Piepenbrink) pointwise inequality; virial identity at the top | the harmonic-oscillator lower bound hλ is built in (sharp in the model); also pointwise certificates | the analyticity domain of U(r(r*)) in complex r* must contain the contours (branch points from r = 0 and r = ∞ images); cross terms from the variable phase; contour design per energy window | yes | 4–6 sessions; risk **medium-high** (contour design) |
| **M4** hyperboloidal energy and commutator identities (P4 S4-2/S4-3 transfer) | redshift and commutation at the horizons (H² level) are a direct analogue of P4's sonic-point treatment, exponents affine in the spectral parameter in both cases | at the trapped set a real multiplier gives only a shrinking strip; a uniform strip needs the commuted unstable-manifold (ladder) structure, which turns into M3/M2 | partially | 5–8 sessions; risk **high** for uniformity; not recommended as primary |
| **M5** cite semiclassical theorems ([SBZ97], [BH08], [WZ11], [Dya12], [Dya16], [Ian17], [Stu26]) | the pseudo-pole lattice, the essential gap γ₀, finiteness of the resonances above −γ₀+ε | **every threshold is "for l large enough"**; making one explicit means redoing its proof (M1 or M3) | no | not a route; supplies the qualitative input of [HV17] (α < γ₀) |

### 4.1 What l₀ to expect (planning heuristics, not results)

If the final relative error at the top is C/(l+½) (M1-quality), exclusion
needs about C/(l₀+½) ≤ (1 − c*)/2, so l₀ ≈ 2C/(1−c*). Olver-type constants
for smooth analytic problems are typically C ~ 2–20. That gives:
- l₀ ≈ 6–60 at P_B (T1);
- l₀ ≈ 9–90 for the two-sided P2-6;
- l₀ ≈ 12–120 at P_A.

If the error is only C/(l+½)^{1/2} (naive local-quadratic models), then
l₀ ≈ (2C/(1−c*))²: about 35–3,500 at P_B and 145–14,000 at P_A. That is
outside the part-(iii) budget (§5). **The plan therefore aims at M1-quality
O(1/l) error control.** M2's centre or M3's contour must be built from the
uniform (Liouville-mapped) approximation.

**Verdict (i): achievable in principle with explicit constants and this
toolset** (Arb, pointwise certificates, Taylor models, closing boxes). It is
new mathematics: an explicit barrier-top theory with complex energy. Estimated
4–8 sessions to a certificate, with a P4-like risk of a residual h → 0 corner.

Subjective odds:
- at P_B: the S2 probe passes, ~55%; a full P2-1 certificate within about
  three quarters, ~40–50%;
- at P_A: ~20–30%.

## 5. Part (iii): sizing

**Method.**
- Frobenius series at r₊ and r_c with exponent −iω/(2κ), as ω-Taylor models
  with tail bounds (P4's sonic-series pattern; bounds in the style of [Mez19],
  re-implemented in python-flint because ore_algebra is Sage-only).
- Taylor-series continuation of (R) with step length a fraction of the distance
  to the nearest singularity. More steps are needed near r₊ at P_B, where
  r₊ − r₋ = 0.097 against 0.194 at P_A; the count grows only logarithmically.
- Jost-normalized W_l.
- Winding by telescoped argument increments on dyadic tiles (P4 `modecount`
  design).
- An analyticity cover (W_l has poles only at E₊ ∪ E_c, far from R_l).
- Krawczyk enclosure of the l = 0 mode at P_B (P2-6).

**Float evidence (S0).**
- **Independent reproduction of [CCDHJ18] Table I** (l = 0 NE; l = 1 dS;
  l = 1, 2, 5, 10 PS) to |Δ| ≤ 4.6e−11, which is the printed precision.
- Float winding numbers of W_l around R_l (Im ∈ [−0.55, 0.1]κ₋):
  **1 for l = 0 and 0 for l = 1, 2, 5, at both points**, with
  min |W/ω| on the contours ≥ 1.25 (P_B) and ≥ 0.41 (P_A).
- The local phase rate |d arg(W/ω)/d Re ω| on the bottom edge is at most
  2.84 rad/κ₋ at P_B and 5.4 rad/κ₋ at P_A, flat in l from 2 to 20. Tile
  widths of 0.5κ₋ (P_B) and 0.3κ₋ (P_A) are therefore adequate.
- The sub-barrier dynamic range is about 10^{1.7(l+½)}, so precision must
  grow to roughly 1.7·l decimal digits.
- The float cost per W evaluation grows about linearly in l (0.09 s at l = 0,
  1.75 s at l = 200).

**Tiles and cost.** The contour per l is half the rectangle, using the symmetry
ω → −ω̄ (the full rectangle for l = 0). Tiles per l are about 24.6(l+½) + 10.6
at P_B and 18.4(l+½) + 17.7 at P_A. Per-tile cost is planned at c₀(1 + l/10)
with c₀ = 3 s (range 1–10 s; P4's heavier 4D tiles took about 40 s). On one
16-vCPU node with 14 workers:

| l₀ | P_B tiles | P_B wall time (range) | P_A tiles | P_A wall time (range) |
|---|---|---|---|---|
| 20 | 5.1 k | 0.7 h (0.2–2.3 h) | 4.0 k | 0.5 h (0.2–1.8 h) |
| 50 | 31 k | 8 h (2.6–27 h) | 24 k | 6 h (2–20 h) |
| 100 | 124 k | 2.3 d (0.8–7.8 d) | 94 k | 1.8 d (0.6–5.9 d) |
| 200 | 494 k | 17 d (6–58 d) | 372 k | 13 d (4–44 d) |

**Parallelization.** The work is embarrassingly parallel over (l, tile). There
is no shared background: the coefficients are exact algebraic balls,
precomputed once per worker in under 1 s. Memory is well under 1 GB per worker,
results are checkpointed per l, and the merge is serial and trivial (P4 ran 290
tiles on 8 workers in 21 min). Speed-up is near-linear to 14–16 workers. **The
budget line is l₀ ≲ 100**, which is days on one node. l₀ ≈ 200 means weeks,
and beyond that the large-l asymptotics of part (i), not part (iii), must carry
the load.

## 6. Risks (ranked)

1. **Part (i) does not close at the needed constants.** The open items are the
   complex-energy barrier top with O(1/l) explicit errors and its h → 0 closing
   corner, which is P4 S4-3's open-sliver pattern. If the projected l₀ exceeds
   about 100, part (iii) becomes unaffordable even when (i) is true.
   Mitigations: the P_B point (c* = 0.33), three independent routes
   (M1/M2/M3), and the S2 kill gate.
2. **Point choice and approval.** Without v1.2, P_A is the only target. Part
   (i) must then reach 0.67 γ₀ against a gap attained only as l → ∞. P_A's
   binding family (photon sphere) is also the numerically fragile one [CDP23].
   Risk (1) roughly doubles.
3. **The conversion chain and the cost tail.**
   - The equivalence lemma D1′ (radial Frobenius QNMs = [Vas13]/[HV17]
     resonances, with multiplicity) and the quantitative hypotheses of
     [HV17] §§2.5–2.6 (0 simple, α < γ₀) must be audited in writing at S1.
     A missed hypothesis would void P2-5 without touching P2-4.
   - Separately, part-(iii) precision and step growth at large l, and
     near-extremal stiffness at P_B, may push c₀ toward 10 s and move the
     table above by 3×.

Further risks:
- Stack gap: the D-finite bounds tool [Mez16] is Sage-only, so the series
  bounds must be re-implemented in python-flint (P4 precedent; about 1
  session).
- Scientific value: the scalar field is the toy model. The physical case is
  gravito-electromagnetic [DRS18], which is out of scope, though the same
  machinery applies later.

## 7. Recommendation and proposed v1.2 revision-log entry

**Recommendation: GO with a modified entry point.** The proposed entry text
below is for your approval. **The .tex has not been edited.**

> **v1.2** (proposed 2026-10-16; P2 S0 gate, `problems/P2/notes/s0-gate.md`).
> **P2 Status.**
> (a) The open item is a resonance-free region
> {Im ω ≥ −(½+δ)κ₋}∖{0} for some δ > 0, with ω = 0 a simple resonance with
> constant resonant states and (½+δ)κ₋ below the trapping gap
> γ₀ = (2r_P)⁻¹(F(r_P)(2 − 3M/r_P))^{1/2}. These are the hypotheses under which
> [hintzvasy2017] gives H^{1+δ−0} regularity at the Cauchy horizon. The open
> strip {−κ₋/2 < Im ω ≤ 0}∖{0} yields only β ≥ ½.
> (b) Conjecture 1.5 of [hintz2025] is the near-extremal-limit form of this
> statement. At fixed parameters the requirement is the explicit strip. At the
> representative point of [cardoso2018, Table I] the binding family is the
> photon-sphere family in the eikonal limit (β = γ₀/κ₋ ≈ 0.764), not the
> near-extremal family.
> (c) [dias2018] concerns Kerr–de Sitter (scalar and gravitational
> perturbations). Coupled gravito-electromagnetic perturbations of RNdS show a
> stronger violation for smooth data [drs2018].
> (d) Bibliographic: [mezzarobba2016] is an arXiv extended abstract of an
> ICMS 2016 talk; the full truncation bounds are M. Mezzarobba, Ann. Henri
> Lebesgue 2, 99 (2019).
>
> **P2 Entry point.**
> - Replace "for an explicit set of subextremal RNdS parameters in the regime
>   where numerics put β comfortably above 1/2" by "at the pre-registered points
>   ΛM² = 3/50 with Q/Q_max = 999/1000 (primary: the near-extremal window of
>   [cardoso2018, Fig. 3], where the dominant mode is an isolated ℓ = 0
>   near-extremal mode with β ≈ 0.91, so that part (i) is needed only on a strip
>   of width ≈ 0.33γ₀) and Q/Q_max = 249/250 (second: [cardoso2018, Table I],
>   β ≈ 0.764, attained only as ℓ → ∞, so part (i) is needed up to ≈ 0.67γ₀)".
> - Replace "(i) an explicit large-ℓ exclusion (barrier-top WKB bounds with
>   explicit error constants)" by "(i) an explicit large-ℓ exclusion, uniform
>   over the whole strip (above-barrier, barrier-top, sub-barrier and
>   low-frequency regimes) for ℓ ≥ ℓ₀ with explicit ℓ₀".
> - Add: "At the primary point the target includes a two-sided certified
>   enclosure of β (an interval-Newton enclosure of the ℓ = 0 near-extremal
>   mode)."
>
> **Bibliography.** Add [drs2018] O. J. C. Dias, H. S. Reall, J. E. Santos,
> "Strong cosmic censorship: taking the rough with the smooth," J. High Energy
> Phys. 10 (2018) 001, arXiv:1808.02895.

If v1.2 is declined, FORMULATION v0 falls back to P_A as the sole target (its
Revision log says so). In that case the gate verdict should be read as GO
**only at reduced odds** (§4).

**Staged plan with kill criteria.**
- **S1** (1–2 sessions). Arb parameters; the D1′ lemma; the P2-0 sympy
  identity; KAT-1…KAT-7, including an independent hyperboloidal spectral code
  for P_B; certified P2-2 thresholds W_l for l ≤ 50.
  *Exit:* all KATs pass.
- **S2** (2–3 sessions): **part-(i) feasibility probe.** Certify the
  barrier-top window at P_B, c* = 0.33, at the discrete values l = 20, 50, 100,
  using the inequalities of the leading part-(i) route (M2 with an
  M1-quality centre; M3 as the alternate) rather than part-(iii) zero counting.
  Fit the error constant C.
  *Kill:* if no route certifies l = 100, or the projected l₀ is above 100 at
  P_B, stop. Record P2 as NO-GO-for-now, with the certified finite-l and
  part-(ii) results logged as partial evidence and explicitly *not* the entry
  point.
- **S3** (2–3 sessions plus hours to days of compute). The part-(iii) pipeline,
  runs for l < l₀ in parallel, and the l = 0 Krawczyk enclosure.
- **S4** (3–5 sessions). Part (i) uniform in h, including the h → 0 closing
  boxes and the low-frequency throat regime.
- **S5**. Assembly, adversarial review, ledger rows P2-0…P2-6, then P2-7 (P_A)
  if part (i) reaches 0.67γ₀.

Total: about 9–15 sessions plus 1–3 days of 16-vCPU compute, if l₀ ≤ 100.

## 8. Reproduction

```
.venv/bin/python problems/P2/src/s0_sanity.py geom      # constants, fractions of gamma0     (<1 s)
.venv/bin/python problems/P2/src/s0_sanity.py modes     # Table I reproduction + P_B spectrum (~20 s)
.venv/bin/python problems/P2/src/s0_sanity.py rate      # phase rates on the bottom edge      (~15 s)
.venv/bin/python problems/P2/src/s0_sanity.py winding   # float zero counts in the boxes R_l  (~95 s)
```

All outputs are floating point. The cost-scaling probe (W evaluation time at
l = 0…200) was a scratch run; its numbers are quoted in §5.
