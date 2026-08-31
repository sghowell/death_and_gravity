# P9(b) S1 — derivation chain, convention pins, known-answer tests

Status: S1 complete, 2026-08-31.  Code: `problems/P9/b/derivation/` (exact
sympy, P4 house style: every step an identity test), tests in
`problems/P9/b/tests/`, numeric KGB check in `problems/P9/b/src/p9b/`.
Run: `uv run pytest problems/P9/b/tests -q` from the repo root.

## 1. What was derived from what

Machinery (`tools.py`, `ibp.py`, `geometry.py`): single scalar Fourier mode
(every perturbation eps*f(t)*cos(kx)), exact Taylor jets to O(eps^2), exact
x-averaging, and an integration-by-parts canonicaliser that returns an
explicit total-derivative certificate F with L = L_canon + dF/dt re-verified
by sympy each time.  ADM geometry (3d Christoffels/Ricci, extrinsic
curvature, full 4d Ricci scalar) is computed from the metric by index loops —
GLV14 (73)-(75) are outputs, not inputs.

Chain (all exact, all tested):
1. `actions.py`: the FROZEN action = GLV14 eq. (86) transcribed term by term
   (five-alpha ADM quadratic action; alpha_B in GLV convention); GLV14
   eq. (87) (GPV operator form) transcribed separately; matter = minimally
   coupled k-essence P(Y) (dust = c_m^2 -> 0 limit).
2. `background.py`: eps^1 (tadpole) variation of [(87) + matter] gives the
   background equations — GLV14 (88)-(89) reproduced EXACTLY — plus the
   matter eom d/dt(a^3 sigma' P_Y) = 0, hence sigma'' = -3 H sigma' c_m^2 and
   rho_m' + 3H(rho_m+p_m) = 0.  Family-B definitions (FORMULATION 1.3) and
   the DE conservation identity.  Unitary-gauge k-essence operators derived:
   c = X P_X, Lambda = X P_X - P, M_2^4 = X^2 P_XX.
3. `dictionary.py` (deliverable v): [(87)]_eps^2 with (88)-(89) substituted
   equals [(86)]_eps^2 under the pinned dictionary (GLV14 table, verified):
       M^2 = M*^2 f + 2 m4^2,           M^2 H^2 aK = 2c + 4 M_2^4,
       M^2 H aB_GLV = (M*^2 fdot - m3^3)/2,   M^2 aT = -2 m4^2,
       M^2 aH = 2(mt4^2 - m4^2),        M^2 H aM = d(M^2)/dt,
   up to an exact d/dt and the matter-coupling remainder
       a^3 [ (3/2) rho_m zeta dN - (9/4) p_m zeta^2 ]   (x-averaged),
   which is the background-equation-proportional piece GLV absorbed when
   deriving (86); it vanishes in vacuum.  f(R) limit (deliverable vi):
   c = M_2^4 = m3^3 = m4^2 = mt4^2 = 0 gives alpha_K = 0 and
   alpha_B^BS = -alpha_M exactly (digest D10 confirmed).
4. `reduction.py` (deliverable ii): the total quadratic action (86)+matter
   (with the remainder above) has dN, psi algebraic; they are eliminated
   exactly.  Vacuum: GLV14 (77) [dN = zetadot/(H(1+aB)), exact for all five
   alphas], (79)-(80) [kinetic M^2 alpha/(1+aB)^2 with alpha = aK + 6 aB^2 =
   (83); gradient = (80) verbatim, k-independent, zero mass term], and BS14
   (3.12) normalisation Q_S = 2 M*^2 D/(2-aB_BS)^2 under the map.  With
   matter: 2-field (zeta, dsigma) UV kinetic/gradient matrices; speeds from
   det(c^2 T + a^2 U/k^2) = 0.
5. `rows.py`: I1 (deliverable iii) and the k-essence dictionary (vii).
6. `kgb.py`: INDEPENDENT route — covariant L = K(X) + G(X) box(phi) in
   unitary gauge, expanded from scratch (box(phi) included).  Background:
   rho = 2X K_X - K - 3 G_X H phid^3, p = K + G_X phid^2 phidd; Bianchi gives
   rhodot + 3H(rho+p) = phid (Jdot + 3HJ), J = phid(K_X - 3 H phid G_X) —
   DPSV (40)-(41) under the signature map G_DPSV = -G_ours.  The eps^2 action
   equals frozen (86) with M^2 = MPl^2, aT = aH = 0 and derived
       H^2 M^2 aK = 2X(K_X + 2X K_XX) - 12 phid X H (G_X + X G_XX),
       H M^2 aB_BS = -2 phid X G_X       (= BS14 (A.7) with G_3 = -G),
   an exact zero-residual match (this simultaneously re-verifies the BS map
   on a concrete family).

## 2. Convention pins verified (equation numbers as fetched 2026-08-31, ar5iv)

- GLV14 (arXiv:1411.3712): (86) five-alpha action [frozen object]; (82)
  alpha_H definition; (83) no-ghost alpha_K + 6 alpha_B^2 > 0; (84)-(85)
  c_s^2; (68) alpha_M = H^-1 dln M^2/dt; (73)-(75) perturbation variables;
  (77) lapse constraint; (79)-(81) reduced action; (87) GPV-form action;
  (88)-(89) background equations.  All reproduced exactly by the chain.
- BS14 (arXiv:1404.3713): (3.4) alpha_M; (3.12)-(3.13) Q_S, D, c_s^2 —
  reproduced exactly under alpha_B^BS = -2 alpha_B^GLV (deliverable iv); the
  map is forced: the same reduction matches both papers' formulas.
- GPV (arXiv:1210.0201): operator names in the frozen dictionary follow
  GLV14 (87)/GLPV 1304.4840.  1210.0201's own names: mbar_1^3 = m3^3,
  Mbar_2^2 = -Mbar_3^2 = 2 m4^2, mu_1^2 = mt4^2/2, f(t) identical.
- DPSV (arXiv:1008.0048): model (59) K = -X, G = mu X in signature (+,-,-,-);
  shift charge (40)-(41).  Signature map to our conventions: G -> -G.

## 3. Discrepancies / clarifications found (none silently adapted)

- GLV14 internal notation: (73) says N^i = delta^{ij} d_j psi but (74) holds
  for the lower-index potential N_i = d_i psi_l with psi_l = a^2 psi.  We use
  (73) literally; psi is auxiliary, nothing downstream changes.  Recorded.
- GLV14 (85) with matter is NOT the exact UV eigen-speed of the coupled
  (zeta, dsigma) system when alpha_H != 0.  Derived exact statements:
    * product rule (all five alphas free, any c_m):
          c_1^2 c_2^2 = c_m^2 * c_s^2(85)     [exact];
    * alpha_H = 0: speeds are exactly {c_m^2, c_s^2(3.13) = c_s^2(85)};
    * dust (c_m^2 -> 0, p_m = 0): one mode exactly soft; the DE eigen-speed is
          c_s^2 = c_s^2(85)|dust + alpha_H^2 rho_m / (M^2 H^2 alpha),
      i.e. matter term -(1+2 alpha_H) rho_m/(M^2 H^2 alpha) instead of
      (85)'s -(1+alpha_H)^2 rho_m/(M^2 H^2 alpha)  (kinetic matter mixing,
      cf. D'Amico et al. arXiv:1609.01272).  FORMULATION D2's "c_s^2 derived
      from the frozen action" is pinned to this exact eigen-speed (Revision
      log v1.1); it agrees with the published pins whenever alpha_H = 0
      (12 of 16 ladder rows).
- FORMULATION 1.2 wrote "m̄4^2" for the operator GLV14 (87) calls m4^2, and
  attributes the (M_2^4, m_3^3, m4^2, mt4^2) names to 1210.0201; naming
  clarified in the Revision log (dictionary content unaffected, verified).

## 4. b0 anchor and dictionaries (verified)

- I1 (S = {alpha_K}, family B, exact):  alpha_K c_s^2 = 3 Omega_DE (1+w_DE)
  = (-2 MPl^2 Hdot - rho_m)/(MPl^2 H^2).  Test: `test_I1_exact`.
- k-essence (deliverable vii): with derived c = X P_X, M_2^4 = X^2 P_XX,
  dictionary alpha_K = (2X P_X + 4X^2 P_XX)/(M^2 H^2), and derived Friedmann
  Hdot = -(2X P_X + rho_m + p_m)/(2 M^2):  I1 <=> c_s^2 = P_X/(P_X + 2X P_XX),
  and rho_DE + p_DE = 2 X P_X.  Vikman's statement in frozen conventions.

## 5. KGB numeric known-answer (digest D2; not a witness)

`src/p9b/kgb_check.py`, DPSV model (mu = 1, MPl^2 = 1, dust), ICs
phid_0 = 0.5, Omega_m0 ~ 0.9, integrated 8 e-folds with the DERIVED
background system and the DERIVED D, c_s^2:
  crossing at N_c ~ 0.194 (w: -0.982 -> -1.013, transversal), rho_DE > 0,
  D >= 0.526 and c_s^2 >= 0.462 on the bracketing window (N in [0, 0.594]);
  w -> -1^- at the J = 0 attractor; Friedmann constraint drift < 1e-6.
Stable phantom crossing with {alpha_K, alpha_B} reproduced as DPSV claim.

## 6. S2 go/no-go: the 8 no-alpha_B rows

The derived exact dust eigen-speed gives, for alpha_B = 0 (D = alpha_K),
exact row targets (test `test_row_targets_no_alphaB`; below f-hat =
3 Omega_DE (1+w_DE) evaluated with the row's M^2):
  {K}      aK c_s^2 = 3 Omega_DE (1+w_DE)                       [pure sign no-go: b0 READY]
  {K,T}    aK c_s^2 = 3 Omega_DE (1+w_DE) - 2 alpha_T
  {K,M}    aK c_s^2 = 2 alpha_M - 2 Hdot/H^2 - rho_m/(M^2 H^2)
  {K,H}    aK c_s^2 = 2 aH + 2 aHdot/H - 2(1+aH) Hdot/H^2 - (1+2 aH) rho_m/(MPl^2 H^2)
  (combined rows: sums of the offsets, same machinery.)
All 8 targets are well-formed polynomial identities in the alphas, their
first derivatives, and background data: S2 is GO.  Honest flags for S2:
only {K} is a pure sign-forced no-go.  {K,T} has a bounded rescue budget
-2 alpha_T <= 2 (c_T^2 >= 0): near the crossing f -> 0, so a small negative
alpha_T can make the target positive through the crossing — this row is an
(E)-candidate, not an obvious (N).  {K,M} and {K,H} have offsets (2 alpha_M;
free aHdot) unbounded by the frozen |alpha| <= 10 bounds (C^1 bounds only) —
genuinely open, as FORMULATION already flags; the S2 no-go attempts there
need extra stated hypotheses or will flip to (E)/CONJECTURED.

## 7. Test inventory

test_derivation.py: geometry pins (74)-(75), ADM inverse, R4 background;
(88)-(89); matter eom + conservation; family-B conservation; k-essence
operators; GPV dictionary + remainder; f(R).   test_reduction.py: (77),
(79)-(80)+(3.12), (83); (3.13) factorisation; (85) product rule; dust mixing
term; I1; k-essence; row targets.   test_kgb.py: box, rho/p, DPSV (40)-(41),
quadratic match, numeric crossing.  All exact except the numeric check.
