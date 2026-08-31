# P9(b) literature digest

Compiled 2026-08-31 (web-checked that day; arXiv/journal ids verified by
search+fetch unless marked [unverified detail]). Cited from FORMULATION.md
as D-numbers. Quotes kept ≤ 5 lines per source.

## 1. The single-field no-go lineage (input to b0)

- **D1. Vikman 2005**, "Can dark energy evolve to the phantom?",
  Phys. Rev. D 71, 023515 (2005), astro-ph/0407107. For single-field
  k-essence P(φ,X) minimally coupled in FLRW: crossing w = −1 is impossible
  for stable perturbations; transitions can only pass through points with
  P_X = 0 (degenerate/unstable), i.e. crossing is "hard to realize only in
  the presence of linear terms in X". Our b0 is this statement re-proven as
  one exact α-basis identity (I1) from the frozen action.
- **D3. Quintom review: Cai–Saridakis–Setare–Xia**, "Quintom cosmology:
  theoretical implications and observations", Phys. Rept. 493 (2010) 1–60,
  arXiv:0909.2776. States and proves a "no-go theorem" for a single perfect
  fluid or single scalar with generic P(φ,X): stable crossing forbidden;
  two-field quintessence+phantom realizations carry an explicit ghost.
  Their proof is perturbation-divergence at crossing; ours is the algebraic
  sign identity — same content, certificate-shaped.
- **D3'. 2025–26 quintom updates (the "Cai et al. 2026 update" exists).**
  (i) Qiu–Cai–Liu–Li–Evslin–Zhang, "A Focused Review of Quintom Cosmology:
  From Quintom Dark Energy to Quintom Bounce", arXiv:2511.19994 (Nov 2025);
  (ii) Ren–Li–Liu–Cai–Li–Zhang, "A short review on Quintom dark energy
  theory", arXiv:2606.01360 (RAA invited, May 2026, rev. Jul 2026) — both
  restate the no-go and survey crossing realizations post-DESI-DR2. Neither
  abstract advertises a minimal-operator-content classification theorem
  [full-text check of both is an S1 to-do].

## 2. Crossing constructions (input to b1 existence rows)

- **D2. Deffayet–Pujolàs–Sawicki–Vikman**, "Imperfect dark energy from
  kinetic gravity braiding", JCAP 10 (2010) 026, arXiv:1008.0048. The
  braiding class L = K(φ,X) + G(φ,X)□φ: the scalar "is able to cross the
  phantom divide with neither ghosts nor gradient instabilities" (their
  abstract). This is the expected (E) winner {α_K, α_B}; known-answer test
  §3.2(3) reproduces one of their crossings.
- **D6. Creminelli–D'Amico–Noreña–Vernizzi**, "The effective theory of
  quintessence: the w < −1 side unveiled", JCAP 02 (2009) 018,
  arXiv:0811.0827. EFT of single-field DE perturbations: "no general
  pathology associated to w < −1 or in crossing the phantom divide";
  stability on the w < −1 side requires behaving as k-essence with
  "virtually zero speed of sound", or the extrinsic-curvature (braiding-type
  m₃³ ↔ α_B) operators. This is the EFT-language ancestor of the b1 dichotomy:
  either c_s² → 0 degeneracy (excluded by our frozen bounds) or braiding on.
- **D7. Matsumoto**, "Phantom crossing dark energy in Horndeski's theory",
  Phys. Rev. D 97, 123538 (2018), arXiv:1712.10015. Constructs Horndeski
  crossing models ("general solution to overcome the difficulty of phantom
  dark energy"); model-building, not a minimality classification.
- **D8. Cataneo–Koyama**, "Non-parametric exploration of minimally coupled
  gravity with phantom crossing", Phys. Rev. D 113, 122006 (2026),
  arXiv:2512.13691. Luminal KGB (mochi_class stable-EFT basis), ghost and
  gradient stability enforced from the outset; finds "viable
  phantom-crossing solutions exist without conformal coupling". Empirical
  scan — explicitly no theorem about minimal operator content (checked).
- **D9. Post-DESI Horndeski crossing constructions (context only):**
  Naidoo–Hallam–Baker–Sirera, "Constraints on Horndeski Gravity with Phantom
  Crossing", arXiv:2606.20794 (asymptotic cubic Galileons, minimally
  coupled); "Cosmology after Phantom Crossing by Horndeski Gravity",
  arXiv:2512.03139; braided DE with momentum exchange, arXiv:2607.26447;
  "Quintessential dark energy crossing the phantom divide", arXiv:2508.19101.
  All are model constructions/fits, none states a classification theorem.

## 3. EFT-of-DE bases and the standard stability conditions (b0/b1 targets)

- **D4. Gubitosi–Piazza–Vernizzi**, "The effective field theory of dark
  energy", JCAP 02 (2013) 032, arXiv:1210.0201. Unitary-gauge action from
  unbroken spatial diffs; background fixed by three operators {f(t)R, Λ(t),
  c(t)g⁰⁰}; perturbation operators M₂⁴(δg⁰⁰)², m₃³ δK δg⁰⁰, m₄²-type
  (δK² − δK^μ_ν δK^ν_μ), etc. Source of the operator dictionary in
  FORMULATION §1.2. Companion: Gleyzes–Langlois–Piazza–Vernizzi, "Essential
  building blocks of dark energy", JCAP 08 (2013) 025, arXiv:1304.4840.
- **D5. Bellini–Sawicki**, "Maximal freedom at minimum cost: linear
  large-scale structure in general modifications of gravity", JCAP 07 (2014)
  050, arXiv:1404.3713. Defines {α_K, α_B, α_M, α_T} (α_M: eq. (3.4);
  others: eqs. (A.7)–(A.10)); quadratic action eq. (3.12); stability
  eq. (3.13):

      Q_s = 2M*²D/(2−α_B)² > 0,   D ≡ α_K + (3/2)α_B²  > 0,
      c_s² = −[(2−α_B)(Ḣ − ½H²α_B(1+α_T) − H²(α_M−α_T)) − Hα̇_B
              + (ρ̃_m+p̃_m)] / (H²D) ≥ 0,   ρ̃ ≡ ρ/M*²,

  (transcribed from the ar5iv fetch 2026-08-31; b0/b1 are proven against the
  re-derivation from the frozen action, with (3.13) as known-answer test —
  any transcription slip dies in S1). Setting α_B = α_T = α_M = 0 and using
  Friedmann gives identity I1: α_K c_s² = 3Ω_DE(1+w_DE) — the b0 engine.
- **D5'. Gleyzes–Langlois–Vernizzi**, "A unifying description of dark
  energy", Int. J. Mod. Phys. D 23 (2015) 1443010, arXiv:1411.3712 [journal
  ref unverified]. Five-α ADM quadratic action, their eq. (86); α_H defined
  eq. (82); no-ghost eq. (83): α_K + 6α_B² > 0 (their α_B convention,
  = −α_B^BS/2, consistent with D of D5); c_s² with α_H terms eq. (85).
  Beyond-Horndeski origin of α_H: Gleyzes–Langlois–Piazza–Vernizzi,
  "Healthy theories beyond Horndeski", PRL 114, 211101 (2015),
  arXiv:1404.6495, and "Exploring gravitational theories beyond Horndeski",
  JCAP 02 (2015) 018, arXiv:1408.1952.
- **D10. Convention cross-checks:** hi_class (arXiv:1605.06102 /
  1909.01828) implements the BS convention; f(R): α_B = −α_M, α_K = 0
  [convention detail to re-verify in S1]. EFTCAMB stability guide:
  Frusciante et al., arXiv:1601.04064.

## 4. Positivity (input to b2) — and the honest gap

- **D11. Adams–Arkani-Hamed–Dubovsky–Nicolis–Rattazzi**, "Causality,
  analyticity and an IR obstruction to UV completion", JHEP 10 (2006) 014,
  hep-th/0602178. Flat-space, Lorentz-invariant, massive 2→2 forward-limit
  positivity. Does NOT directly apply on FLRW: the DE background
  spontaneously breaks boosts, the scalar is effectively massless, and no
  S-matrix exists on the expanding background. This gap is the reason b2 is
  conditional.
- **D12. Melville–Noller**, "Positivity in the sky", Phys. Rev. D 101,
  021502(R) (2020), arXiv:1904.05874. Tree-level flat-space bounds for
  (shift-symmetric) Horndeski, eqs. (16)–(17), e.g.
  2Ḡ₂,X Ḡ₄,X ≥ −Ḡ₂,XX Ḡ₄ and 2Ḡ₄,XX + 2Ḡ₄,X²/Ḡ₄ ≤ Ḡ₃,X²; example-model
  α-form eq. (10): α_B ≤ 2α_T/(1+α_T). Their own caveat: bounds derived for
  massive particles on flat backgrounds, applied on FLRW by *assuming* they
  "continue to hold (at least approximately)". This assumption, stated
  verbatim, is our frozen POS assumption set.
- **D13. Melville–Noller**, "Positivity bounds on dark energy: when matter
  matters", JCAP (2022), arXiv:2103.06855 [journal ref unverified]: matter
  loops shift the bounds; relevant if POS is ever upgraded.
- **D14. Boost-breaking positivity (closer to the right framework, not yet
  a usable DE bound):** Grall–Melville, "Positivity bounds without boosts",
  Phys. Rev. D 105, L121301 (2022), arXiv:2102.05683;
  Creminelli–Janssen–Senatore, "Positivity bounds on effective field
  theories with spontaneously broken Lorentz invariance", JHEP 09 (2022)
  201, arXiv:2207.14224 (assumes conformal UV, uses retarded Green's
  functions — no S-matrix needed). de Rham–Melville–Tolley lineage:
  "Improved positivity bounds and massive gravity", JHEP 04 (2018) 083,
  arXiv:1710.09611 [id from memory, verify in S1]; survey: de Rham et al.,
  "Snowmass White Paper: UV constraints on IR physics", arXiv:2203.06805.
- **D15. dS/cosmological positivity 2024–26:** "Propagator positivity
  bounds for cosmological correlators", arXiv:2512.20706, PRD (2026) —
  two-sided bounds for *heavy* fields on dS (wrong regime for a light DE
  scalar, but the first genuinely-dS tower); Horndeski/DHOST triple-crossing
  bounds, arXiv:2306.06639; "Phenomenology of Horndeski gravity under
  positivity bounds", arXiv:2403.13096; multiple-vacua bounds,
  arXiv:2202.01222; EFT-of-perturbations positivity, arXiv:1908.08644.
  Verdict for b2: no theorem yet transports any of these to a bound on
  (α_K, α_B, α_T, α_M, α_H)(t) on our background at our regularity — hence
  POS-conditional formulation, with D14/D15 as the upgrade path.

## 5. Is (b1) already in the literature? (prior-art check, plain answer)

Partially, as folklore and examples; not as a classification theorem.
- Necessity direction: D1/D3 prove {α_K} fails (our b0 = certificate-grade
  restatement). D6 shows in EFT language that stability at/near crossing
  needs braiding-type operators or c_s² ≈ 0. No source states or proves
  "the minimal sets S ⊆ {α_K,…,α_H} admitting stable crossing are exactly
  {…}" over an explicit ladder with explicit boundedness hypotheses.
- Existence direction: D2 (KGB) is a stable-crossing existence proof by
  example (analytic + numeric, not validated numerics); D8/D9 give
  post-DESI stable-EFT scans and constructions.
- Rows genuinely undecided in the literature: α_M-without-α_B,
  α_H-without-α_B (and their combinations) — no stable-crossing example and
  no no-go found in the searches above.
- Our delta: (i) b0 as an exact machine-checked identity in frozen
  conventions; (ii) the 16-row ladder decided with certificates on both
  sides (Arb-validated witnesses; exact no-go identities), including the
  currently undecided rows; (iii) the positivity intersection as an explicit
  conditional semialgebraic statement rather than a prior in an MCMC (D12's
  usage). If S1's full-text pass over D3', D8, D9 finds a genuine
  classification claim, this section is revised and the delta restated.

## 6. DESI DR2 status (motivation only; no data enters (b))

- **D16. DESI collaboration**, "DESI DR2 results II: measurements of baryon
  acoustic oscillations and cosmological constraints", arXiv:2503.14738.
  w₀w_aCDM preferred over ΛCDM at 3.1σ (DESI+CMB), 2.8–4.2σ with SN
  compilations; favored quadrant w₀ > −1, w_a < 0; best fit crosses w = −1
  at z ≈ 0.4–0.5 [crossing-z from secondary summaries, not pinned]. Dataset-
  consistency disputes are live (e.g. arXiv:2506.15091 "Could we be fooled
  about phantom crossing?", arXiv:2506.19053); the theory question stands
  regardless (pre-registration §P9).

## 7. Open items carried into S1

1. Full-text pass: D3' (both reviews), D8, D9 — confirm no hidden
   classification theorem; pin the DPSV example used in test §3.2(3).
2. Re-verify every equation pin (§1.2 of FORMULATION) against the papers,
   including the α_B convention map and MN eq. (10) sign convention.
3. Verify the [unverified detail] items above: 1710.09611 id, journal refs
   for 1411.3712 and 2103.06855, f(R) α_B = −α_M convention.
4. Decide the POS variant (D12 alone vs. D12+D13) before freezing b2.
