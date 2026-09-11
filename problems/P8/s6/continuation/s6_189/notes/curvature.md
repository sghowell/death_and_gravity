# Direct dynamic and spatial tensor curvature

Use the unchanged clock chart with N=1,shift0,u=t and
g_ij=-a^2(exp gamma)_ij. The spatial gamma is symmetric,
transverse and traceless. Its volume is exactly a^3,
X=R_target=1,trace K=3H,uHu=0 and full affine source
S=0. R_target denotes the parent coefficient, not
the curvature scalar. The fixed scalar profile and
original homogeneous M1 are unchanged and have no
gamma dependence on this chart.

The independent four-dimensional calculation uses

    g_mu_nu=diag(1,-a^2 exp q,-a^2 exp(-q),-a^2),
    q=q(t,z),

and contracts the full connection and dynamic Riemann
tensor before the quadratic expansion. It gives exactly

    R_P8=-6(H'+2H^2)
                  -(q_t^2-a^-2 q_z^2)/2.

Thus for a general compact TT field the quadratic
scalar-curvature integral against a time-only
coefficient is one quarter of the Frobenius kinetic
minus gradient form, modulo a spatial divergence.
Its first TT variation is zero. The plus matrix used
above has Frobenius norm squared2.

The same direct calculation gives the quadratic Weyl
scalar for that plus component:

    C^2_2=(q_tt+Hq_t+a^-2 q_zz)^2-4a^-2 q_tz^2.

Let Dq=q_tt+Hq_t-a^-2 q_zz. The literal density
difference is the compact boundary

    a^3[C^2_2-(Dq)^2]
      =4{partial_t(a q_t q_zz)-partial_z(a q_t q_tz)}.

This proves the coefficient, including the mixed
spatial/time derivatives, without inferring it from
the vanishing background Weyl tensor.

For all spatial momenta, the homogeneous isotropic
background makes the quadratic form diagonal in
Fourier momentum. Rotate each momentum to z. The
cross polarization follows by a constant transverse
rotation, and a reflection makes the mixed plus/cross
term vanish. Both orthonormal TT polarizations have
the same coefficient. Independent rotated linearized
four-index curvature tests check both simultaneously.
Plancherel and compact integration by parts then give

    integral sqrt(-g) C^2 at quadratic order
       =1/2 integral a^3 tr[(D gamma)^2],
    D=partial_t^2+H partial_t-a^-2 Delta.

Equivalently this follows from four-dimensional
conformal invariance and the flat compact Euler
identity. The direct curved calculation and its
boundary supply independent verification here.
Zero spatial momentum requires no inverse momentum
operation in the resulting local operator.
