# Fresh scalar reference and the finite prescription

Use FK signature (+---), R=6(H'+2H^2), R_00=3(H'+H^2), G_00=-3H^2.
Define the two local tensors by their complete covariant expressions

    H3_ab=R_a^c R_cb-(2/3)R R_ab
             -(1/2)g_ab R_cd R^cd+(1/4)g_ab R^2,
    I_ab=2R R_ab-(1/2)g_ab R^2
             +2 nabla_a nabla_b R-2g_ab Box R.

The scalar conformal reference is fixed by

    T_conf=hbar/(2880pi^2)[H3+beta_S I].                  (1)

beta_S is a NEW finite real R^2-prescription coordinate. It is not
identified with a photon coefficient or any earlier scalar gamma.
The Lambda/Newton choices are separate; finite Weyl variation
vanishes in this conformally flat background. No independent
nonzero gravitational curvature-squared term is silently discarded:
if present its contribution must enter the additional-source budget.

The scalar Euler input is a4=1/360 and the conformally flat reference
stress is Eq. (23) of
[Herzog and Huang](https://arxiv.org/pdf/1301.5002).
Their metric and Ricci tensor convention map simultaneously to
opposites of FK, making H3_HH,00=-H3_FK,00. Their factor
-2a4/(16pi^2) therefore gives +H3_FK/(2880pi^2).
The zero-type-D choice is beta_S=0, not a universal prescription.

For clarity, the complete actual FLRW contractions, with the common
hbar/(2880pi^2) removed, are derived here:

    rho=3H^4+beta_S(18H'^2-108H^2 H'-36H H''),
    trace=12H^2(H^2+H')
          -36beta_S[H'''+4H'^2+7H H''+12H^2 H'],
    pressure=(rho-trace)/3,
    E=-3(H^4+2H^2 H')
         +beta_S[18H'''+90H'^2+90H H''+108H^2 H'].

The arbitrary-jet identities in stress.py check the Ricci contractions,
the I variation, its trace -6 Box R, physical conservation
rho'+3H(rho+pressure)=0, and E=rho-trace/2. This fixes the state
reference using (1), not by integrating conservation with an
unfixed a^-4 radiation constant.

For H=1/(2t), I_E vanishes and E_conf=hbar/(5120pi^2 t^4).
In the beta_S=0 prescription rho_conf=hbar H^4/(960pi^2).
These are scalar values, distinct from the photon values.

Given nonnegative jet caps v=(v0,v1,v2,v3), the triangle inequality
bounds the absolute E numerator by

    V_S(v)=3(v0^4+2v0^2 v1)
       +|beta_S|(18v3+90v1^2+90v0 v2+108v0^2 v1).

Each monomial has weight four under proper-clock rescaling. This
is the only reference-loss polynomial used in the new cosmological
cost. The trace anomaly is retained even though target-minus-reference
stress is traceless.
