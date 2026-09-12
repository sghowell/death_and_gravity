# Exact synchronous maps without singular inversions

The physical ADM chart is g00=N^2-h(beta,beta), g0i=-(h beta)i, gij=-hij, with N=1+n and h=a^2 exp(Q). Direct evaluation of the covariant metric Lie derivative for xi=(eta,chi) gives

    n_xi=eta',
    beta_xi=chi'-a^-2 grad eta,
    Q_xi=2H eta I+grad chi+(grad chi)^T.

For source G define eta_G=Iminus n_G and chi_G=Iminus(beta_G+a^-2 grad eta_G), where Iminus integrates from the unchanged left endpoint. Then G-Gxi has zero lapse/shift and spatial component Qsyn_G=Q_G-2H eta_G I-symgrad chi_G. Its gauge vector vanishes on the same initial neighborhood. A smooth extension of the flow outside the observation slab is identity near preparation; the pulled-back state has exactly the same initial covariance.

For the detector use Iplus f=-integral_t^right f, whose derivative is again f. This gauge vector vanishes on the final neighborhood but need not vanish initially. The two choices enter different Ward identities.

Neither map divides by H, momentum, bounce energy or a constraint symbol. Here a>=1. They are regular at t=0 and P=0 and cost at most two spatial derivatives of lapse, one of shift and time primitives. They are not a gauge-fixed dynamical inverse.
