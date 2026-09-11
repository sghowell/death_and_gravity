# Covariant chart, including the null-gradient surface

Write R for the unchanged positive S6.109 tensor coefficient, not the
Ricci scalar. Let X=g_physical^(mu nu)u_mu u_nu, with P8 signature+---.
That source proves1/2<R<6/5 on every real u and -1/4096<X<6/5 and
R-1=X^2 A(u,X), with A analytic. Its scalar and curvature clock jets
are independently re-read here, not replaced by a global clock polynomial.

Define the algebraic change of metric variables

    g_physical=C g_hat+D du du,
    C=R^(-1/2), D=(1-C)/X.

Here R is evaluated at u and X_hat. Since C+DX=1, inversion gives
X_physical=X_hat and

    g_hat=C^(-1)g_physical+(1-C^(-1))du du/X.

The determinant ratio is C^3=R^(-3/2). The ten-dimensional Jacobian
at fixed du is a rank-one update of C times identity. Its remaining
eigenvalue is C-X C_X-X^2 D_X=1; its determinant is C^9, never zero.
These are explicit inverses, not only a unitary-gauge assertion.

Both apparently divided coefficients are analytic through X=0:
1-C=(R-1)/[sqrt(R)(sqrt(R)+1)]. In particular D=O(X), including
nonzero null gradients. Lorentz signature is preserved: on timelike
or spacelike gradients the orthogonal three directions are rescaled
positively and the gradient norm is unchanged; the null case follows
continuously, with nonzero determinant. This does not identify the
hat cone with the physical matter cone. Matter remains coupled to
g_physical; there is no quotient by a matter-changing transformation.

For X=s^2>0 in scalar unitary gauge, the lapse is unchanged and
h_physical=R^(-1/2)h_hat. Set omega=-log(R)/4. The physical trace is

    K_phys=K_hat+3s omega_u+6s omega_X V, V=n.partial s.

The exact target coefficients are
B=-R/2, C_ADM=-s R_X, D_ADM=-3X R_X^2/(4R) and
E_ADM=-X R_X+7X^2 R_X^2/(4R).
Direct substitution removes K_hat V and V^2. The remaining linear
coefficient is L_V=3s^2 R_u R_X/(2R), with volume R^(-3/4).
Take I_s=R^(-3/4)L_V and I(u,1)=0. Its full time/spatial integration
by parts contributes -I K_hat-s I_u. Dropping this boundary would
change the lapse constraint.

The spatial-curvature conformal transformation and its spatial
boundary give exactly E_ADM+2B4 delta_c^2+
4(B4-2X B4_X)delta_c=0, where B4=R/2 and delta_c=X R_X/(2R).
Thus neither lapse velocities nor spatial lapse derivatives remain
in this light-action chart. These are generic-R identities, not
specializations to the old linear-in-X coefficient.
