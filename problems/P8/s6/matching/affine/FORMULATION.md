# S6.37: a source-preserving auxiliary metric-affine lift of CD/M1

This is a new classical auxiliary-parent investigation, not a modification
of any earlier action or certificate. The target is the original
[CD/M1 witness](../../../certificates/witness-CD_matter.json), with its
original free chi, physical metric and all scalar coefficient functions.
The adopted [S6 V/G/B contract](../../FORMULATION.md) is unchanged.

## Domain and conventions

Use dimensionless M=tau=1 variables initially. The entire closed target
tube is u real and X_repo in [9/10,11/10], within an open neighborhood.
The source-paper convention is (-+++): g_source=-g_repo, x=-X_repo,
R_source=-R_repo. This is a signature dictionary, not a reassignment of
the matter metric. All lower connection indices remain independent;
the derivative index in Gamma^a_bc is the last one. Only the four
projective directions Gamma^a_bc -> Gamma^a_bc+delta^a_b U_c are gauge.

Let h=(1+u^2)^3, w=h-1-x, and choose

    f=w/(2h), p=sqrt(w/h)/2,
    c=2(p-f)/x, J4=(1/2-f)/x^2.

The parent action in source conventions is

    integral sqrt|g| [p R_Gamma+c G_Gamma^{ab} phi_a phi_b
                     +J2+J3 L3_Gamma+J4 L4_Gamma - (d chi)^2/2].

Here G_Gamma is the literal double-dual affine Einstein contraction;
L3_Gamma and L4_Gamma are the two-epsilon projectively invariant
Galileon terms, not versions obtained by prematurely integrating by
parts or raising a differentiated scalar index. Their definitions and
the curvature boundary are specified in [the action proof](notes/connection.md).

## Exact matching statement

Define q(u,x) as the unique solution of

    q_x+[1/(2x)-3p_x/(2p)]q=3(p_x/p) f_u,  q(u,-1)=0.

It is a smooth local coefficient function, not a history-dependent
field inverse. A positive integrating factor gives a definite integral
over x inside the same tube. Choose

    J3=-(q+f_u)/(4px),
    J2=F_repo(u,-x)+x q_u-3x(q+2f_u)^2/(16p^2).

The denominator of J3 is the product 4*p*x, not a derivative. F_repo is
the full frozen scalar function, not only its clock value or X jets.

On this tube the full affine Euler equation has one solution modulo
projective gauge. Eliminating the connection yields precisely the
original CD/M1 action up to an explicit covariant divergence. In
source conventions its coefficients are

    f=w/(2h), alpha1=alpha2=0, alpha3=1/(hx),
    alpha4=(11x-4h+4)/(4hxw), alpha5=-1/(hxw).

The lower-order terms satisfy Q1=q, Q2=2q_x and P-xq_u=F_repo(u,-x).
The boundary identity is

    q Box(phi)+2q_x phi^a phi_ab phi^b
      =div(q grad(phi))-xq_u.

All five second-derivative operators, the curvature coefficient, the
full scalar function and the free matter action are matched on an open
neighborhood of the closed tube. A solution of the original target
remaining inside this domain therefore lifts, locally and modulo the
projective gauge, to a solution of the auxiliary parent. This follows
by varying the exact reduced action, not by imposing only a background
or tensor equation. In particular the original CD/M1 rolling solution
has the same physical metric and matter backreaction.

## Full algebraic rank and a quantitative bound

For a timelike scalar rest frame with x=-s^2, keep all 64 distortion
components. Trace gauge kappa^a_ac=0 leaves 60 components. The literal
quotient Hessian determinant in the declared coordinates is

    -2^52 p^72(8p^2-1)^3.

On the closed tube, p^2 is in [9/40,11/40], so 8p^2-1>=4/5. Its
row-sum inverse norm is at most 880109/36000<25. This is a coordinate
norm in an orthonormal scalar rest frame and the displayed trace gauge,
not a Lorentz-invariant propagating mass. The compact source denominator
Delta=2p-cx+2J4*x^2=1 alone would miss an additional degeneracy at
p^2=1/8: the quotient rank there is 57 despite compatible scalar forcing.
That locus is not in the actual CD tube. It is not promoted to a kinetic
ghost diagnosis or a universal obstruction to other extensions.

The paper's printed Q2 bracket p-3p_x is a retained negative control.
Direct elimination of the literal action gives p-3x*p_x. A separate
pure p(phi,x) R_Gamma calculation reproduces the same correction.
No silent repair of the source formula is an input to the connection solve.

## Physical units and verification boundary

For M^2,tau>0, multiply p,c,J4 and the effective quadratic coefficients
by M^2; J3 and q by M^2/tau; J2 by M^2/tau^2. The action in rescaled
coordinates has prefactor M^2*tau^2. Thus Delta_physical=M^2, not one.
The original physical clock is phi=tau*u and free chi=M*chi_bar.

The written covariance, algebraic elimination and smooth-ODE proofs are
checked by exact symbolic identities, a separately authored coefficient
and signature audit, full connection Euler tests and rejected-domain
controls. This is not proof-assistant formalization.

No connection kinetic term, heavy gap, controlled kinetic promotion,
vacuum extension, loop bound, Regge/IR remainder, positivity verdict or
UV completion is established. Exact auxiliary equivalence is not itself
completion of B, V or G. Scoped earlier P8(a)/(b) results are unchanged;
original P8 remains open.
