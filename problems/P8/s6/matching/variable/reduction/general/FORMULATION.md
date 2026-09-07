# S6.31: regular-root universality and genuine next-order freedom

This is a new formal-action theorem, not a change to the frozen S6.20
example or its background. It extends the stationary-metric calculation
of S6.30 to general ordinary interaction coefficients and distinguishes
a genuine change of action from a change of operator representative.

## General class and retained-order statement

Keep the physical metric g, signature +---, curvature R_B, canonical
clock phi and canonical free matter chi. The Einstein densities have
constant positive squared Planck masses MG2 and MF2, respectively.
The interaction is -2 sqrt(g) sum beta_n(phi) e_n(sqrt(g^-1 f)), with
all five beta_n smooth on an open clock interval. No extra curvature,
derivative interaction, disformal matter coupling or metric redefinition
is part of this class.

Assume a smooth positive proportional root r(phi) of the hidden metric's
own algebraic equation on the interval. Define

    F_root = beta1+3r beta2+3r^2 beta3+r^3 beta4 = 0,
    P = beta1+2r beta2+r^2 beta3 != 0.

The partial derivative of F_root with respect to r on the root is
-3P/r. Thus this is the simple-root domain, not an assumed vacuum or
solution of the physical-g and clock equations. P can have either sign;
invertibility is not a health or mass-sign assertion. A nonsimple root
P=0 has vanishing full quadratic potential Hessian in this chart and
requires a different expansion, not division by a small number set to zero.

Writing f=r^2(g+h), H=g^-1 h, the full symmetric-metric quadratic
potential is -rP([H^2]-[H]^2)/4. Its ten-component Hessian determinant
in the stated Lorentz chart is 3(rP)^10/16. The own-f solution through
two formal derivatives and its action through four are

    h2 = (MF2*r/P)(B-g tr_g(B)/3),  B=G[r^2 g],
    L4 = kappa (B:B-tr_g(B)^2/3),
    kappa = MF2^2*r^3/(4P).

Here MF2 is already a squared mass. The zero- and two-derivative action
is the explicit variable-conformal Einstein action with Planck
coefficient MG2+MF2*r^2, clock kinetic coefficient 1/2-3MF2*r_phi^2,
unchanged free chi, and potential
-2(beta0+3r beta1+3r^2 beta2+r^3 beta3).

After the full variable-coefficient integrations by parts in S6.30,
the specified scalar-tensor normal form still obeys A1=2F2_X and
A3=A4=A5=0. The separate variable curvature-square operator remains.
Consequently this retained-order action cannot exactly match the
original CD coefficient tube in the same physical metric by a
scalar-only clock change. This extends only the retained-order
statement, not an all-order obstruction or a controlled error estimate.

## A genuine deformation invisible through retained order

At fixed r and P, any smooth eta(phi) defines a change of the parent
interaction coefficients

    Delta(beta0,beta1,beta2,beta3,beta4)
      = (3r eta,-2eta,eta/r,0,-eta/r^3).

It preserves F_root, P and the stationary potential polynomial. Thus
it preserves every formal stationary-action coefficient through degree
four, the hidden Einstein Hessian, and h2. It does not claim to preserve
the S6.20 rolling solution of the changed full action.

The cubic potential depends also on Z=beta1+3r beta2+2r^2 beta3:

    V3 = r{Z[H]^3-3(P+Z)[H][H^2]+(3P+2Z)[H^3]}/24.

The deformation changes Z by eta. By stationarity, f4 terms cancel
from the degree-six action. The complete formal degree-six action
difference between the two parents is therefore exactly

    Delta L6 = r eta ([h2]^3-3[h2][h2^2]+2[h2^3])/24.

The Einstein second-variation term is unchanged. The full degree-six
action itself is not computed here. A compact conformal variation on
a constant-curvature off-shell metric proves that this action
difference is not generically a boundary. This control is not an
assertion that the chosen test metric solves either parent.

## Potential table and representative control

The full curvature-linear part of V3[h2] is given as an exact covariant
table, with all clock gradients and Hessians retained. It is explicitly
potential-only. A declared projection onto its Hessian-free curvature
operators is not an invariant after the higher-operator complement is
changed by integration by parts.

In particular, div(W(phi) R_B X grad(phi)) contains both W_phi R_B X^2
and three retained curvature-derivative/Hessian operators. A coefficient
W that vanishes at the center can have a nonzero derivative there.
The resulting named F2_X shift is canceled by the full complement in
the action's Euler equation. A literal tensor boundary test checks this
cancellation and rejects the center-first projection shortcut.

## Explicit nonclaims

No functional inverse, homogeneous state, rolling spectral gap,
convergent expansion, all-operator or loop bound, healthy global
background, positivity/V/G test or full classification is supplied.
Neither the Z deformation nor an apparent projected Xi supplies a
controlled match to CD. Other nonminimal/derivative parents are not
covered by this theorem. Original P8 remains open.
