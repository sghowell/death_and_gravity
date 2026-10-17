# P2 literature digest — strong cosmic censorship with Λ > 0 (RNdS spectral gap)

Status: S0, 2026-10-16. Every arXiv id below was fetched on 2026-10-16 (arXiv
abstract page or source tarball). Journal references were taken from the arXiv
record, from INSPIRE (`inspirehep.net/api`), or from the bibliography of a paper
read in full. In the last case the entry says "(per X bbl)". Entries marked
**unverified** were not checked independently. Equation, figure, table and
theorem numbers are the cited paper's own, read from its TeX source unless the
entry says otherwise. Abbreviations in brackets are used in `FORMULATION.md` and
`s0-gate.md`.

Section 8 re-checks the P2 paragraph of
`docs/problems/open-problems-theoretical-cosmology-2026.tex` (v1.1) and lists
the deltas.

---

## 1. The claim: β > 1/2 for near-extremal RNdS (non-rigorous numerics)

**[CCDHJ18] V. Cardoso, J. L. Costa, K. Destounis, P. Hintz, A. Jansen,
"Quasinormal modes and Strong Cosmic Censorship", Phys. Rev. Lett. 120, 031103
(2018), arXiv:1711.10502** (v2, 2018-01-11; TeX source read in full).
- Model: a massless minimally coupled scalar on fixed RNdS,
  F = 1 − 2M/r + Q²/r² − Λr²/3 (eq. 3); κ_* = ½|F′(r_*)| (eq. 4). The radial
  equation is d²Φ/dr*² + (ω² − V_l)Φ = 0 (eq. 6) with
  V_l = F(l(l+1)/r² + F′/r) (eq. 7). QNM boundary conditions
  Φ ~ e^{−iωr*} at r₊ and Φ ~ e^{+iωr*} at r_c (eq. 8), time dependence e^{−iωt}.
- Definitions (eq. 9): α ≡ inf_{l,n}{−Im ω_{ln}}, β ≡ α/κ₋. Footnote 2: "For
  l=0 there is a zero mode, corresponding to Φ₀ in Eq. (1), which we ignore
  here." Eq. (2): the field is H^{1/2+β} at the Cauchy horizon, citing [HV17]
  and [CF17].
- Three families. Photon-sphere (PS) modes: the dominant one is "determined by
  the large l limit (and n=0)", and "l=10 or l=100 provide good approximations
  of the imaginary parts". de Sitter (dS) modes deform the pure-dS modes
  ω/κ_c^{dS} = −il (n=0) and −i(l+n+1) (n≥1) (eqs. 11–12), with
  κ_c^{dS} = √(Λ/3). Near-extremal (NE) modes approach
  ω_NE = −i(l+n+1)κ₋ = −i(l+n+1)κ₊ in the extremal limit (eq. 13). [Hin25]
  shows eq. (13) is wrong for l ≥ 1.
- **Table I** (supplementary material): the dominant ω₀/κ₋ for
  l = 0, 1, 2, 3, 4, 5, 10, 100 at ΛM² = 0.06, Q/Q_max = 0.996, from Jansen's
  spectral code at (N, p) = (400, 200) and (450, 225). The text quotes 61–82
  agreed digits but prints 10 decimals. Values used as known answers:
  l=0 −0.8539013779 i (NE); l=1 3.2426164126 − 0.7958326323 i (PS) and
  −1.5003853731 i (dS); l=2 5.4796067815 − 0.7754179185 i;
  l=3 7.7016152057 − 0.7699348317 i; l=4 9.9181960834 − 0.7676996545 i;
  l=5 12.1322536641 − 0.7665731858 i; l=10 23.1897597770 − 0.7649242108 i;
  l=100 222.0602900249 − 0.7643094446 i. The caption says "the true dominant
  mode being that with l → ∞", so at this point β is attained only in the
  eikonal limit.
- **Fig. 3** (`Figure_nearExtremal.pdf`, rendered and read): four panels at
  ΛM² = 0, 0.02, 0.06, 0.14. The tick labels give β = 1/2 at
  Q/Q_max = 0.9903, 0.9917, 0.9920, 0.9945, and the NE mode becomes dominant at
  0.9959, 0.9966, 0.9967, 0.9975 respectively. In the ΛM² = 0.06 panel the
  l=0 NE curve stays between about −0.8κ₋ and −1.0κ₋ for Q/Q_max ∈ [0.9967, 1).
  The PS curve plotted there is l = 10, not the l → ∞ limit.
- Fig. 2 is a contour plot of β over the whole parameter space. The dS family
  dominates for ΛM² ≲ 0.02 and the PS family elsewhere, with a thin
  NE-dominated sliver at maximal charge. Fig. 4 shows eigenfunctions at
  ΛM² = 0.02, Q/Q_max ≈ 0.9966: NE −0.87005i, dS −0.87043i and
  PS(l=10) 26.448 − 0.89101i, all in units of κ₋.
- **S0 reproduction**: `src/s0_sanity.py modes` reproduces all six Table I
  entries it tests (l = 0 NE, 1 dS, 1, 2, 5, 10 PS) to |Δ| ≤ 4.6e−11. That is
  the rounding of the printed 10 decimals. The method is independent:
  Frobenius series at r₊ and r_c, Taylor-series continuation of the D-finite
  radial ODE, and Wronskian zeros. This pins the point as exactly
  ΛM² = 0.06, Q = 0.996·Q_max(Λ). The eikonal gap γ₀/κ₋ = 0.7643026660
  ([HV17] §2.7 formula) matches the l = 100 entry, 0.7643094446, up to an
  O(l⁻²) difference of 6.8e−6.

**[DERS18] O. J. C. Dias, F. C. Eperon, H. S. Reall, J. E. Santos, "Strong
cosmic censorship in de Sitter space", Phys. Rev. D 97, 104060 (2018),
arXiv:1801.09694.** This paper is about **Kerr–de Sitter**, with scalar and
linearized gravitational perturbations: "for any non-extremal value of the black
hole parameters, there are quasinormal modes which decay sufficiently slowly to
ensure that strong cosmic censorship is respected." It does not study
gravito-electromagnetic perturbations of RNdS; that is [DRS18].

**[DRS18] O. J. C. Dias, H. S. Reall, J. E. Santos, "Strong cosmic censorship:
taking the rough with the smooth", JHEP 10 (2018) 001, arXiv:1808.02895** (TeX
source read, §§1, 5.3). This is the follow-up on rough versus smooth data, and it
treats coupled gravitational and electromagnetic perturbations of RNdS. It uses
the Kodama–Ishibashi master variables Φ_±, vector and scalar type, which are
isospectral (§4.2).
- It finds the same three families. The near-horizon AdS₂×S² analysis (§5.3)
  gives ω/κ₋ = −i(n + ½ + √(¼ + η̂_±)) with
  η̂_± = 1 + Ξl(l+1) ± √([1 + Ξl(l+1)]² − Ξ²(l+2)(l+1)l(l−1)). The slowest
  NE mode is Φ₋ with l = 2.
- Result: the violation is "much worse". For a sufficiently large near-extremal
  black hole, perturbations from smooth data extend through the Cauchy horizon
  in an arbitrarily smooth way (C^r for any r).
- On the old Brady–Moss–Myers argument for SCC, they show it holds only for
  initial data that is not C¹ at the event horizon, which supports the
  rough-data formulation of [DSR18].
- What this adds for P2: a second, gravitationally realistic field model, which
  is out of scope here (FORMULATION §7). Its master equations have rational
  potentials, so the same D-finite machinery would apply.

**[DRS19] O. J. C. Dias, H. S. Reall, J. E. Santos, "Strong cosmic censorship
for charged de Sitter black holes with a charged scalar field", Class. Quantum
Grav. 36, 045005 (2019), arXiv:1808.04832.** A charged scalar "improves the
situation but there is always a neighbourhood of extremality in which strong
cosmic censorship is violated".

**[CCDHJ18b] V. Cardoso, J. L. Costa, K. Destounis, P. Hintz, A. Jansen,
"Strong cosmic censorship in charged black-hole spacetimes: still subtle",
Phys. Rev. D 98, 104007 (2018), arXiv:1808.03631.** For charged massive scalars
"there still exists a finite region in parameter space where Strong Cosmic
Censorship is expected to be violated".

**[LZCCN19] R. Luna, M. Zilhão, V. Cardoso, J. L. Costa, J. Natário, "Strong
cosmic censorship: the nonlinear story", Phys. Rev. D 99, 064014 (2019),
arXiv:1810.00886.** Nonlinear, spherically symmetric Einstein–Maxwell–scalar
evolutions with Λ > 0 find that mass inflation "indeed does not occur in the
near extremal regime".

**[HWZ20] S. Hollands, R. M. Wald, J. Zahn, "Quantum instability of the Cauchy
horizon in Reissner–Nordström–deSitter spacetime", Class. Quantum Grav. 37,
115009 (2020), arXiv:1912.06047.** The quantum stress tensor diverges as
T_VV ~ C/V² at the Cauchy horizon, so quantum effects restore SCC. This is
outside P2's classical scope.

**[DDG24] A. Davey, O. J. C. Dias, D. Sola Gil, "Strong Cosmic Censorship in
Kerr-Newman-de Sitter", JHEP 07 (2024) 113, arXiv:2404.03724.** A scan of the
KNdS parameter space that locates the boundary between respecting and violating
SCC, for linear scalar fields. It is cited by [Hin25] for matched-asymptotics NE
formulas.

**[CDP23] A. Courty, K. Destounis, P. Pani, "Spectral (in)stability of
quasinormal modes and strong cosmic censorship", Phys. Rev. D 108, 104027
(2023), arXiv:2307.11155.** A small bump in the potential destabilizes PS modes
but leaves the dominant dS and NE modes stable, and the RNdS violation survives.
Relevance to P2: β is a property of the exact operator, but PS-dominated points
are numerically fragile. A further reason to make the NE-dominated point
primary.

**[CCDP25] A. Chrysostomou, A. S. Cornell, A. Deandrea, S. C. Park, "A note on
Strong Cosmic Censorship and its violation in Reissner–Nordström de Sitter
black hole space-times", Class. Quantum Grav. 42, 107001 (2025),
arXiv:2501.12968.** A review together with numerical parameter maps. Not
rigorous.

**[SZ26] Z. Stuchlík, A. Zhidenko, "Nonoscillatory gravitational quasinormal
modes of Reissner–Nordström–de Sitter spacetime", Phys. Rev. D 114, 084038
(2026), arXiv:2506.09829.** Analytic expressions for purely imaginary
gravitational modes. No rigorous claims.

**[AD26] V. Arrigoni, G. Dusson, "Computing resonances of perturbed Schrödinger
equations: Application to Reissner-Norsdtröm-de Sitter black holes"** (title
misspelling as on arXiv), **arXiv:2606.18770** (2026-06-17; TeX source read).
- Method: zeros of Jost-function Wronskians, found with a deflated Newton
  iteration. The reference Jost solutions use Hatsuda's Heun representation
  [Hat20].
- Their potential (their eq. `eq:pot_reissner`) contains an extra +(2Λ/3)F.
  That is the conformally coupled field, a different model from [CCDHJ18],
  and the source of the Heun structure.
- They give no error bounds, and their SCC check covers only l ∈ {0,…,25}.
  It is the nearest numerical analogue of our part (iii), but it is not
  certified.

**[Hat20] Y. Hatsuda, "Quasinormal modes of Kerr-de Sitter black holes via the
Heun function", Class. Quantum Grav. 38, 025015 (INSPIRE lists 2020),
arXiv:2006.08957.**

**[Jan17] A. Jansen, "Overdamped modes in Schwarzschild-de Sitter and a
Mathematica package for the numerical computation of quasinormal modes", Eur.
Phys. J. Plus 132, 546 (2017), arXiv:1709.09178.** QNMspectral: hyperboloidal
(ingoing Eddington–Finkelstein) Chebyshev collocation, solved as a generalized
eigenproblem. This is the code behind [CCDHJ18].

**[JMA21] J. L. Jaramillo, R. Panosso Macedo, L. Al Sheikh, "Pseudospectrum and
black hole quasinormal mode (in)stability", Phys. Rev. X 11, 031003 (2021),
arXiv:2004.06434.** QNM spectra are strongly non-normal in the hyperboloidal
energy norm, and the pseudospectra open up at high frequency. The consequence
for P2: exclusion arguments based on resolvent norms or numerical range in
energy norms are pessimistic. The S0 plan avoids norms in part (iii), where
zeros of an analytic Wronskian are counted.

## 2. Rigorous spectral theory of the exterior (asymptotic, non-explicit constants)

**[SBZ97] A. Sá Barreto, M. Zworski, "Distribution of resonances for spherical
black holes", Math. Res. Lett. 4, 103–121 (1997)** (per [Hin25] bbl; not
re-fetched). Complex scaling for Schwarzschild(-de Sitter). High-energy
resonances are close to the lattice ±(l+½)Ω_P − i(n+½)λ_P. The approximation is
asymptotic, with no explicit thresholds.

**[BH08] J.-F. Bony, D. Häfner, "Decay and non-decay of the local energy for the
wave equation on the de Sitter–Schwarzschild metric", Commun. Math. Phys. 282,
697–719 (2008), arXiv:0706.0350** (journal ref per [Hin25] bbl). Resonance
expansion on SdS and local-energy decay. The l-uniform estimates are
semiclassical, with non-explicit constants.

**[MSBV14] R. Melrose, A. Sá Barreto, A. Vasy, "Asymptotics of solutions of the
wave equation on de Sitter–Schwarzschild space", arXiv:0811.2229** (journal
reference **unverified**).

**[Vas13] A. Vasy (appendix by S. Dyatlov), "Microlocal analysis of
asymptotically hyperbolic and Kerr–de Sitter spaces", Invent. Math. 194,
381–513 (2013), arXiv:1012.4391** (journal ref per [Hin25] bbl). Fredholm
framework on spaces that extend across the horizons, so that QNMs become poles
of a meromorphic family. Radial-point thresholds make QNMs H^s-eigenvalues for
Im ω > −(s−½)κ at each horizon.

**[Dya11a] S. Dyatlov, "Quasi-normal modes and exponential energy decay for the
Kerr-de Sitter black hole", Commun. Math. Phys. 306, 119–163 (2011),
arXiv:1003.6128.** **[Dya11b] S. Dyatlov, "Exponential energy decay for
Kerr–de Sitter black holes beyond event horizons", Math. Res. Lett. 18,
1023–1035 (2011), arXiv:1010.5201.** **[Dya12] S. Dyatlov, "Asymptotic
distribution of quasi-normal modes for Kerr-de Sitter black holes", Ann. Henri
Poincaré 13, 1101–1166 (2012), arXiv:1101.1260.**
- These give a Bohr–Sommerfeld description of the high-frequency QNMs in any
  strip and a resonance-free strip below the real axis for large |Re ω|, valid
  for small a.
- All thresholds are of the form "for |Re ω| large enough", with no explicit
  values.

**[Dya15] S. Dyatlov, "Resonance projectors and asymptotics for r-normally
hyperbolic trapped sets", J. Amer. Math. Soc. 28, 311–381 (2015),
arXiv:1301.5633.** **[Dya16] S. Dyatlov, "Spectral gaps for normally hyperbolic
trapping", Ann. Inst. Fourier 66, 55–82 (2016), arXiv:1403.6401.** The essential
gap equals ν_min/2. There are finitely many resonances above Im = −ν_min/2 + ε
and infinitely many below; constants are not explicit.

**[WZ11] J. Wunsch, M. Zworski, "Resolvent estimates for normally hyperbolic
trapped sets", Ann. Henri Poincaré 12, 1349–1385 (2011), arXiv:1003.4640.**
**[NZ15] S. Nonnenmacher, M. Zworski, "Decay of correlations for normally
hyperbolic trapping", Invent. Math. 200, 345–438 (2015), arXiv:1302.4483.** These
supply the high-energy resolvent bounds in a strip below the real axis that [HV17]
needs, as a qualitative input.

**[HV18] P. Hintz, A. Vasy, "The global non-linear stability of the Kerr–de
Sitter family of black holes", Acta Math. 220, 1–206 (2018), arXiv:1606.04014.**
**[Hin18] P. Hintz, "Non-linear stability of the Kerr–Newman–de Sitter family of
charged black holes", Ann. PDE 4, 11 (2018), arXiv:1612.04489.** Nonlinear
stability for small a and Q. Neither proves anything quantitative about the
spectral gap at near-extremal charge.

**[Ian17] A. Iantchenko, "Quasi-normal modes for de Sitter–Reissner–Nordström
black holes", Math. Res. Lett. 24, 83–117 (2017), arXiv:1407.3654.** Barrier-top
resonances for the massless charged Dirac field on RNdS, generalizing [SBZ97].
Asymptotic expansions with explicit coefficients, but remainders without
constants.

**[Bes20] N. Besset, "Decay of the local energy for the charged Klein–Gordon
equation in the exterior De Sitter–Reissner–Nordström spacetime", Ann. Henri
Poincaré 21, 2433–2484 (2020)** (per [Hin25] bbl; arXiv id not checked). Complete
resonance expansion on subextremal RNdS. Discreteness of the QNM set.

**[HX22] P. Hintz, Y. Xie, "Quasinormal modes of small Schwarzschild–de Sitter
black holes", J. Math. Phys. 63, 011509 (2022), arXiv:2105.02347.** In the
zero-mass limit, shallow QNMs converge to the de Sitter QNMs together with the
Schwarzschild resonances. Rigorous, but a limit statement.

**[Hin-KdS] P. Hintz, "Mode stability and shallow quasinormal modes of Kerr–de
Sitter black holes away from extremality", arXiv:2112.14431** ("Accepted for
publication in J. Eur. Math. Soc." per the [Hin25] bbl). Again a limit
(Λ → 0, or small mass) characterization.

**[Hin25] P. Hintz, "Quasinormal modes of near-extremal
Reissner–Nordström–de Sitter spacetimes", arXiv:2504.01734** (2025-04-02; TeX
source read; no journal reference on arXiv or INSPIRE as of 2026-10-16).
- Notation: r_C < r_e < r_c (Cauchy, event, cosmological), κ_C the Cauchy
  surface gravity (our κ₋), ϰ = (r_c² + 2r_e r_c − 3r_e²)/(r_c² + 2r_e r_c + 3r_e²)
  ∈ (0,1), and λ_l⁺(μ) = ½(1 + √(1 + 4(l(l+1) + r_e²μ)/ϰ)).
- **Thm 1.1**: fix r_e < r_c. As r_C ↗ r_e, the QNMs with |ς| ≤ C₀κ_C,
  divided by κ_C, converge in Hausdorff distance and with multiplicity to
  {0} ∪ {−i(λ_l⁺(0)+n)} (μ = 0). Each limit point has multiplicity 2l+1, and
  these QNMs are purely imaginary for small r_e − r_C.
- He states explicitly that [CCDHJ18] eq. (13) "is inaccurate even to leading
  order" for l ≥ 1.
- Hyperboloidal time t* = t − T(r), T′ = T̃/F (§1.1, label `EqItstar`). As
  printed, T̃ = 2(r−r_e)/(r_c−r_e), which lacks the −1 needed for his stated
  normalization T̃(r_e) = −1, T̃(r_c) = +1; we use the normalized version. A
  QNM is a σ with a smooth resonant state u, that is
  □_g(e^{−iσt*}u) = 0 with u ∈ C^∞.
- **Conj. 1.4** (polyhomogeneous dependence). **Conj. 1.5** ("Shallow QNMs"):
  the same convergence holds for all QNMs with Im ς > −C₀κ_C. This is a
  statement about the **limit r_C ↗ r_e** and would give β → 1.
- Remark 1.3: QNMs of extremal RNdS with negative imaginary part need Gevrey
  or complex-scaling tools ([GW21], [SBZ97], [HZ25], [Stu24]).
- At our points ϰ = 0.8485 (P_A) and 0.8588 (P_B), giving
  λ₁⁺ ≈ 2.11, λ₂⁺ ≈ 3.2. The S0 float NE modes at P_B are −0.909, −1.912 and
  −2.891 (κ₋ units, l = 0, 1, 2), consistent with these limits at finite
  distance from extremality.

**[Joy22] J. Joykutty, "Existence of zero-damped quasinormal frequencies for
nearly extremal black holes", Ann. Henri Poincaré 23, 4343–4390 (2022)** (per
[Hin25] bbl). Conformally coupled mass only.

**[Mav23] G. Mavrogiannis, "Morawetz estimates without relative degeneration
and exponential decay on Schwarzschild–de Sitter spacetimes", Ann. Henri
Poincaré 24, 3113 (2023), arXiv:2111.09494.** Exponential decay by
physical-space multipliers. The rate is not tied to the spectral gap, so this is
not a route to β > 1/2.

**[Stu24] T. Stucker, "Quasinormal modes for the Kerr black hole",
arXiv:2407.04612.** **[Stu26] T. Stucker, "Semiclassical estimates near
threshold energies and resonance counting on Schwarzschild black holes",
arXiv:2607.15145** (2026-07-16). Complex scaling combined with microlocal
analysis. [Stu26] proves a Weyl law in sectors and the "absence of high angular
momentum QNM from a disc whose radius grows linearly with the angular momentum".
Constants are not explicit (abstract-level check).

**[HZ25] M. Hitrik, M. Zworski, "Overdamped QNM for Schwarzschild black holes",
Ann. PDE 11, 28 (2025), arXiv:2406.15924.**

## 3. Conversion theorems: spectral gap ⇒ regularity at the Cauchy horizon

**[HV17] P. Hintz, A. Vasy, "Analysis of linear waves near the Cauchy horizon of
cosmological black holes", J. Math. Phys. 58, 081509 (2017), arXiv:1512.08004**
(TeX source read).
- **Thm 1.1** (non-degenerate RNdS with Q ≠ 0): u = u₀ + u′ with
  |u′| ≤ Ce^{−αt₀}, and u′ ∈ e^{−αt₀}H^{1/2+α/κ₁−0}. Here κ₁ is the Cauchy
  horizon surface gravity, our κ₋.
- The quantitative form is the §2.5 proposition on partial asymptotics
  (`PropRNdSPartialAsymp`) together with the §2.6 theorem on conormal
  regularity (`ThmRNdSPartialAsympConormal`). It takes α > 0 such that 0 is
  the only resonance in Im σ ≥ −α, with α < γ₀ so that the high-energy
  estimates of their §2.4 proposition `PropRNdSSemiclassical2` hold.
- Lemma `LemmaRNdSResonanceSupp` (§2.4): 0 is a simple resonance with constant
  resonant states. Its proof relies on [SBZ97, BH08, Dya11a].
- §2.7 computes the essential gap as
  γ₀ = μ(r_P)ν_min/4 = (2r_P)^{−1}(μ(r_P)(2 − 3M/r_P))^{1/2}, with μ = F and
  r_P the photon sphere (r_P² − 3Mr_P + 2Q² = 0, independent of Λ).
- S0 check: (½+δ)κ₋ < γ₀ at both frozen points (FORMULATION §1.3).
- For P2 this is **the linear conversion theorem**: a certified strip
  {Im ω ≥ −(½+δ)κ₋}∖{0} with 0 simple and (½+δ)κ₋ < γ₀ gives u′ ∈ H^{1+δ−0},
  hence H¹_loc, at the Cauchy horizon for every solution with smooth data.

**[CF17] J. L. Costa, A. T. Franzen, "Bounded energy waves on the black hole
interior of Reissner–Nordström–de Sitter", Ann. Henri Poincaré 18, 3371–3398
(2017), arXiv:1607.01018** (TeX source read).
- Hypothesis: the Price law (Assumption 2.2), meaning exponential decay with
  rate p along H⁺ of the energy flux of every T^lΩ^Iφ.
- Main theorem (Thm 3.1): 2 min{p, κ₊} > κ₋ implies bounded energy up to the
  Cauchy horizon and an extension in C⁰ ∩ H¹_loc.
- The rate p follows from the spectral gap by the resonance expansion.
  κ₋ < 2κ₊ holds at both frozen points (κ₊/κ₋ = 0.677 at P_A, 0.823 at P_B).

**[CGNS15a/b, CGNS17] J. L. Costa, P. M. Girão, J. Natário, J. D. Silva, "On the
global uniqueness for the Einstein–Maxwell-scalar field system with a
cosmological constant", Part 1, Class. Quantum Grav. 32, 015017 (2015),
arXiv:1406.7245; Part 2, Commun. Math. Phys. 339, 903 (2015), arXiv:1406.7253;
Part 3, Ann. PDE 3, 8 (2017), arXiv:1406.7261.** Nonlinear, spherically
symmetric, with exact RNdS data on the event horizon.

**[CGNS18] J. L. Costa, P. M. Girão, J. Natário, J. D. Silva, "On the occurrence
of mass inflation for the Einstein–Maxwell-scalar field system with a
cosmological constant and an exponential Price law", Commun. Math. Phys. 361,
289 (2018), arXiv:1707.08975** (TeX source read, Thm 1.1).
- Hypothesis: a two-sided exponential Price law ∂_vφ(0,v) ≈ e^{−sk₊v} on the
  event horizon, with ρ = k₋/k₊ > 1.
- Mass inflation occurs if s < min{ρ, 2}.
- There is no mass inflation, and an extension with C⁰ metric, L²_loc
  Christoffel symbols and H¹_loc field, if ρ < 9/7 and s > (14/9)ρ.
  [CCDHJ18] footnote 1 translates this to "β > 7/9".
- At P_B: ρ = 1.215 < 9/7 and the float β ≈ 0.909 > 7/9. At P_A: ρ = 1.478 and
  β ≈ 0.764, which fail both conditions.
- The theorem is conditional on the nonlinear Price law.

**[Ros25] F. Rossetti, "Strong cosmic censorship for the spherically symmetric
Einstein–Maxwell–charged–Klein–Gordon system with positive Λ: stability of the
Cauchy horizon and H¹ extensions", Ann. Henri Poincaré 26, 675–753 (2025),
arXiv:2309.14420** (TeX source read, §1).
- Hypothesis: a Price-law upper bound |φ| + |∂_vφ| ≤ Ce^{−sv}, with
  v = (t + r*)/2 near H⁺.
- C⁰ extension always holds.
- No mass inflation, plus an extension with L²_loc Christoffel symbols and
  H¹_loc field, if s > K₋ (equivalently β > ½) and ρ = K₋/K₊ < 2. Both frozen
  points satisfy ρ < 2.
- Conditional on the nonlinear Price law.

**[DSR18] M. Dafermos, Y. Shlapentokh-Rothman, "Rough initial data and the
strength of the blue-shift instability on cosmological black holes with Λ > 0",
Class. Quantum Grav. 35, 195010 (2018), arXiv:1805.08764.** For generic data in
a slightly rough class, the local energy blows up at the Cauchy horizon for all
subextremal RNdS and KNdS parameters. This is the Christodoulou statement in the
rough category. P2 is about smooth data only.

## 4. Frameworks in which QNMs are eigenvalues (for the equivalence lemma)

**[War15] C. M. Warnick, "On quasinormal modes of asymptotically anti-de Sitter
black holes", Commun. Math. Phys. 333, 959 (2015), arXiv:1306.5760.**
Hyperboloidal slices: QNMs are eigenvalues of the generator on H^k, valid for
Im ω > −(k−½)κ. This is the model for our equivalence lemma and for the
H^s-threshold bookkeeping in `s0-gate.md` §2.

**[GW21] D. Gajic, C. M. Warnick, "Quasinormal modes in extremal
Reissner–Nordström spacetimes", Commun. Math. Phys. 385, 1395–1498 (2021),
arXiv:1910.08479.** **[GW24] D. Gajic, C. M. Warnick, "Quasinormal modes on Kerr
spacetimes", arXiv:2407.04098.** Gevrey-regular spaces, needed when κ = 0. P2's
points are subextremal (κ > 0), so the analytic (Frobenius) characterization
suffices.

**[GZ21] J. Galkowski, M. Zworski, "Outgoing solutions via Gevrey-2
properties", Ann. PDE 7 (2021), DOI 10.1007/s40818-021-00094-2,
arXiv:2004.07868.** The Gevrey-2 definition agrees with complex scaling.

## 5. Rigorous computation and explicit asymptotics (tools for parts (i)–(iii))

**[Mez16] M. Mezzarobba, "Rigorous multiple-precision evaluation of D-finite
functions in SageMath", arXiv:1607.01967** (2016-07-07). The comments field
says "5th International Congress on Mathematical Software (ICMS 2016)". The
author's own bibliography lists it as an arXiv technical report, the "extended
abstract of a talk at the 5th International Congress on Mathematical Software"
(HAL hal-01342769). It is the `ore_algebra` numerical-evaluation module:
regular-singular connection problems with rigorous error bounds. It is
**SageMath-only, and SageMath is not in the project's locked stack.**

**[Mez19] M. Mezzarobba, "Truncation bounds for differentially finite series",
Ann. Henri Lebesgue 2, 99–148 (2019), DOI 10.5802/ahl.17.** Tail bounds for
series solutions, covering generalized series at regular singular points. It is
the algorithmic reference for re-implementing the bounds in python-flint, as P4
already does for its sonic and centre series (`problems/P4/notes/s2-validated-expansions.md`).

**[Olv74] F. W. J. Olver, *Asymptotics and Special Functions* (Academic Press,
1974; A K Peters reprint 1997).** Ch. 6 covers Liouville–Green with
error-control functions, for real and complex variables (progressive paths).
Ch. 11 covers simple turning points (Airy) with explicit error bounds.

**[Olv75] F. W. J. Olver, "Second-order linear differential equations with two
turning points", Phil. Trans. R. Soc. Lond. A 278, 137–174 (1975)** (citation
confirmed by web search; paper not re-read). Parabolic-cylinder approximations
uniform through the coalescence of two turning points, with explicit error
bounds, for a real coalescence parameter.

**[Dun21] T. M. Dunster, "Uniform asymptotic expansions for solutions of the
parabolic cylinder and Weber equations", J. Classical Anal. 17 (2021) (JCA
17-06), arXiv:2010.05428.** Parabolic cylinder functions with a large parameter
and complex argument, with explicit error bounds.

**[Dun25] T. M. Dunster, "Asymptotic expansions for solutions of differential
equations having coalescing turning points, with an application to Legendre
functions", arXiv:2504.19405** (final version 2025-11-01; to appear in Stud.
Appl. Math. per the arXiv comments; TeX source read). It treats
w″ = {u²f(a,z) + g(z)}w with u large and **real**, z in a complex domain, and a
**real** parameter a ∈ [a₀, a₁] whose two real simple zeros coalesce at a₀.
Expansions are in parabolic cylinder functions with explicit error bounds. For
P2 the coalescence parameter (the energy ω²/(l+½)²) is **complex**. That case is
not covered by [Olv75] or [Dun25], and closing it is the work identified in
`s0-gate.md` §4 (route M1).

**[Joh17] F. Johansson, "Arb: efficient arbitrary-precision midpoint-radius
interval arithmetic", IEEE Trans. Comput. 66, 1281 (2017).** The project uses
Arb through python-flint 0.9.0.

**[BZ16] B. Barker, K. Zumbrun, "Numerical proof of stability of viscous shock
profiles", Math. Models Methods Appl. Sci. 26, 2451 (2016), arXiv:1601.00837**
(journal ref per `problems/P4/notes/theorem-c-plan.md` §3). The precedent for the
pipeline shape: analytic high-frequency bounds, then a validated winding number
on the compact remainder.

**[GJ26] K. Goldstein, V. Jejjala, "A computer assisted existence proof for a
non-Schwarzschild black hole in Einstein–Weyl gravity", arXiv:2609.21797**
(2026-09-18). Taylor models with directed rounding, convergent horizon series
and Poincaré–Miranda matching. It is the only 2026 computer-assisted proof about
a black hole that the search found, and it treats existence only (no QNMs).

## 6. Search log for 2025–26 rigorous QNM work (2026-10-16)

Queries:
- "rigorous enclosure quasinormal modes interval arithmetic computer-assisted
  proof black hole 2025"
- "quasinormal validated numerics / interval arithmetic / computer-assisted
  spectral gap de Sitter black hole rigorous bound 2026"
- "Hintz conjecture shallow quasinormal modes near-extremal RNdS proof 2025
  2026"
- "explicit resonance-free strip Schwarzschild-de Sitter … explicit constants
  large angular momentum barrier top"
- "strong cosmic censorship RNdS rigorous proof β > 1/2 2025 OR 2026"

Hits:
- Non-rigorous numerics: [AD26] (closest in spirit), [SZ26], [CCDP25], and the
  recurrence-based high-accuracy QNM method arXiv:2503.17325 (authors not
  recorded).
- Rigorous but non-explicit: [Stu26].
- Computer-assisted, but existence only (no QNMs): [GJ26].

Nothing found:
- No rigorous enclosure of any black-hole QNM.
- No explicit-constant resonance-free strip for any de Sitter black hole.
- No proof or disproof of [Hin25] Conj. 1.5.
- No journal publication of [Hin25].

A dedicated search may still have missed a preprint posted in the last weeks
before 2026-10-16.

## 7. What transfers from P4 (in-repo)

- Validated argument-principle counting with complex-parameter Taylor models
  (`problems/P4/notes/s2-theorem-b.md` §3.4: 290 tiles, about 40 s per tile,
  1231 s wall on 8 workers).
- Krawczyk zero enclosures (§3.6) and the analyticity cover (§3.5).
- Series at singular points with exponents affine in the spectral parameter
  (P4's sonic point, where σ(κ) = −σ₀ − (1−σ₀)κ). Its counterpart here is the
  horizon exponent −iω/(2κ₊).
- Energy and commutator certificates with pointwise sign conditions in Arb
  (`s4-energy.md`). The once-commuted certificate is in
  `git show p4/theorem-c-s4-3:problems/P4/notes/s4-commutator.md`. Its open
  item, the τ > 224 near-sonic sliver caused by the coupled scaling ξτ (ξ the
  distance to the sonic point), is the closest in-repo analogue of the h → 0
  barrier-top corner of P2 part (i).

## 8. Re-check of the P2 paragraph of the problems document (v1.1) — deltas

| # | .tex claim (Status/Interface/Entry point) | Finding (2026-10-16) |
|---|---|---|
| 1 | numerics find β>1/2 for near-extremal RNdS [cardoso2018] | Confirmed: Fig. 2, Fig. 3 crossings, Table I. Reproduced at S0 to 1e−10 (section 1) |
| 2 | Kerr–de Sitter appears safe [dias2018] | Confirmed: [DERS18] is KdS, with scalar and gravitational perturbations. **Omission:** gravito-EM perturbations of RNdS [DRS18] give a stronger violation and are not cited |
| 3 | charged and nonlinear cases are contested | Softer than the literature. Charged scalars [CCDHJ18b, DRS19] and nonlinear numerics [LZCCN19] all find that a violation persists near extremality. Only quantum effects [HWZ20] restore SCC |
| 4 | interior settled conditionally: rate βκ₋ ⇒ H^{1/2+β−0} [hintzvasy2017] | Confirmed, with precise hypotheses: 0 simple, no other resonance in Im σ ≥ −α, and α < γ₀ (high-energy estimates) |
| 5 | nonlinear and charged analogues [cgns2018, rossetti2025] | Confirmed. Hypotheses: [CGNS18] needs ρ < 9/7 and β > 7/9 under a two-sided Price law. [Ros25] needs ρ < 2 and β > ½ under a Price-law upper bound. Both are conditional on a nonlinear Price law |
| 6 | rough data: Christodoulou holds for all subextremal [dsr2018] | Confirmed (CQG 35, 195010) |
| 7 | Hintz located the NE family, correcting [cardoso2018] for l ≥ 1 | Confirmed (Thm 1.1; "inaccurate even to leading order" for l ≥ 1) |
| 8 | "β>1/2 … reduces to … Conjecture 1.5 of [hintz2025]" | **Imprecise.** Conj. 1.5 is a near-extremal-limit statement (fixed r_e < r_c, r_C ↗ r_e) and would give β → 1. At a fixed point the requirement is the explicit strip. At the [CCDHJ18] Table I point the binding family is PS in the eikonal limit: γ₀/κ₋ = 0.7643 by the [HV17] §2.7 formula, not NE |
| 9 | open item: strip {−κ₋/2 < Im ω ≤ 0}∖{0}, uniform in l | **Imprecise.** The open strip only gives β ≥ ½. The requirement is {Im ω ≥ −(½+δ)κ₋}∖{0} for some δ > 0, with 0 simple and (½+δ)κ₋ < γ₀. The closed upper half-plane is classical (an exact energy identity) |
| 10 | "the radial equation is Fuchsian, so [mezzarobba2016] applies directly" | True: regular singular points at r_n, r₋, r₊, r_c, ∞ (five; Heun-class only for the conformally coupled field). The implementation is SageMath, outside the stack. Re-implementation in python-flint follows the P4 series pattern |
| 11 | bib: [mezzarobba2016] "in Mathematical Software – ICMS 2016 (Springer, 2016)" | **Not confirmed.** The record is an extended abstract of an ICMS 2016 talk, listed by its author as a technical report. The full bounds paper is [Mez19] |
| 12 | bib: [hintz2025] arXiv:2504.01734 (2025) | Correct, still arXiv-only |
| 13 | Entry (i): "barrier-top WKB bounds with explicit error constants" | **Under-specified.** The exclusion for l ≥ l₀ must be uniform over the whole strip: above-barrier, barrier-top, sub-barrier, and low-frequency (near-horizon or NE) regimes |
| 14 | Entry: "Parts (i)–(ii) … do not yet exist in the literature" | Confirmed by the search in section 6. Part (ii) is close to existing technology ([Olv74] Ch. 6) |
