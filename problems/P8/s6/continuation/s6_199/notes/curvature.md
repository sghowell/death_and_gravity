# Four-dimensional curvature and tensor weights

For a minimal scalar, modulo its separately tracked total derivative,

    a4_scalar=R^2/72+(Riem^2-Ric^2)/180.

For the minimal one-form Laplace operator, the rank-four scalar term
is supplemented by tr(RE/6), tr(E^2/2), and tr(Omega^2/12), where the
endomorphism is -Ric and the connection curvature is the one-form
Riemann action. These give respectively -R^2/6, Ric^2/2, -Riem^2/12:

    a4_vector=4a4_scalar-R^2/6+Ric^2/2-Riem^2/12.

The Proca vector-minus-scalar determinant therefore has

    a4_Proca=-Riem^2/15+29Ric^2/60-R^2/8
            =13C^2/120-7Euler/40+R^2/72,

with C^2=Riem^2-2Ric^2+R^2/3 and
Euler=Riem^2-4Ric^2+R^2. This algebra retains Euler. It does not
continue the identity away from four dimensions or discard evanescent
contacts when recovering a finite action.

The vector-minus-scalar construction and its difference from massless
Maxwell theory are explicit in section 4.3, equations (4.24)-(4.25),
of [Gorbar and Shapiro, 2003](https://arxiv.org/pdf/hep-ph/0303124).
This is supporting context; the full stress amplitudes and retarded
normalization here are derived independently.

For a transverse Fourier metric direction represented in COM by the
symmetric spatial matrix H, direct linearized curvature contraction gives

    Riem_quad=s^2 tr(H^2),
    Ric_quad=s^2 (tr(H)^2+tr(H^2))/4,
    R_quad=s^2 tr(H)^2,
    C_quad=s^2 (tr(H^2)-tr(H)^2/3)/2.

Mixed variation in two independent directions yields the Hessians
s^2 P2 and 6s^2 P0. Three independent fixtures reconstruct the full
four-dimensional Riemann tensor and contract it with the Lorentz metric,
rather than substituting the displayed scalar formulas.

It follows that

    32 pi^2 lim_(s->infinity) rho2(s)/s^2 =13/120,
    32 pi^2 lim_(s->infinity) rho0(s)/s^2 =6/72.

The scalar channel retains the massive longitudinal ultraviolet
contribution. These are four-dimensional tensor-weight comparisons,
not a determination of finite local subtraction constants.

Section V of [Franchino-Vinas et al., 2018](https://arxiv.org/pdf/1812.00460)
provides an independent external-metric Proca form-factor context.
Its conclusions specify a second-order curvature expansion on
asymptotically flat space. That expansion is not a proved remainder
estimate for the full curved CD response.

Both external-metric normalization factors give 4/kappa. The spin0
projector is not a canonically reduced propagating scalar and this
calculation supplies no full mixed constraint-reduced norm. S6.193's
finite dimensionally continued prescription remains exactly unchanged.
