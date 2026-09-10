# Interaction-preserving canonical limit

The displayed action uses dimensionless coordinates x/tau
and kappa=M^2 tau^2. Set Phi_bar=sqrt(kappa)u,
Y_bar=kappa X. Physical Phi=Phi_bar/tau. Hold tau fixed:
the light mass is 1/tau, the derivative quartic coefficient
is lambda tau^4 and the DHOST quartic coefficient is
-2gamma tau^6. These are fixed canonical interactions.

For fixed n, sending kappa to infinity would remove the
DHOST term. That valid but different limit is not used for
the retained-interaction matching. Instead take even n to
infinity with kappa_n=n/gamma, fixed nonzero gamma and lambda.
Write eps=kappa^-1/2, u=eps Phi_bar, X=eps^2 Y_bar.
Then on fixed compact canonical field/jet sets,

    S_n=exp[(gamma/eps^2)log(1-eps^4 Y_bar^2)]
       =1-gamma eps^2 Y_bar^2
          +(gamma^2/2)eps^4 Y_bar^4+O(eps^6),
    R-1=-gamma eps^2 Y_bar^2+O(eps^4),
    R_X=-2gamma Y_bar+O(eps^2).

The stated positive domain contains such compact sets with
Y_bar>-1/(4gamma). Its lower X endpoint shrinks in original
units; no fixed open original-X neighborhood for every n is
being silently assumed. Uniform finite derivatives on compact
canonical sets follow from the analytic log expansion and
its convergent local derivative series.

Independent limits of the full lower action and DHOST
coefficients give

    L_limit=Y_bar/2-Phi_bar^2/2
            +lambda Y_bar^2-gamma Phi_bar^4/3
            -2gamma(L3_bar-L4_bar).

In particular the actual induced quartic terms before
co-scaling are lambda Y_bar^2-n Phi_bar^4/(3kappa)
and -2n(L3_bar-L4_bar)/kappa. The switch and localizer terms
are retained in the independent full-function expansion.
There is no cubic scalar interaction at the vacuum.
The potential Phi_bar^2/2+gamma Phi_bar^4/3 is positive.

A canonical graviton has dimensionless curvature O(eps).
The nonconstant tensor factor is O(eps^2); hence its action
mixing is kappa*O(eps^2)*O(eps)=O(eps). Covariant-derivative
graviton insertions in the surviving scalar quartic are also
suppressed by eps. The Einstein self-interactions decouple
in the usual canonical expansion. The unchanged physical
free chi and ordinary Proca sectors do not acquire a direct
flat-vacuum interaction with Phi in this limit. This is a
classical canonical action/finite-vertex limit, not a
nonperturbative interchange of a loop integral with M->infinity.

With v_mu=d_mu Phi, H_mu nu=d_mu d_nu Phi,

    div[X(v^mu box Phi-H^(mu nu)v_nu)]
      =2(L3-L4)+X[(box Phi)^2-H_mu nu H^(mu nu)].

Third derivatives cancel by their commuting flat indices.
The remaining identity is independently checked as a
polynomial in arbitrary Lorentzian first/second jets.
Thus the retained interaction is the flat quartic Galileon
up to the displayed divergence, not a quadratic higher-time-
derivative pole at the vacuum. None of this certifies
nonlinear health on the off-clock transition or finite-M
quantum positivity.
