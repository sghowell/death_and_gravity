# S6.30: source-preserving formal reduction and the CD coefficient defect

This is a new bounded S6 matching audit of the unchanged S6.20 parent.
It does not alter that action, the original CD/M1 witness, or the adopted
S6 V/G/B contract. The result concerns the formal stationary-f action
through four total derivatives, not a proved low-energy expansion.

## Action, branch and physical quantities

Keep both positive Einstein coefficients M^2, signature +--- and the
curvature convention R_B=-6(Hdot+2H^2). The interaction is
-2 sqrt(g) sum beta_n(phi) e_n(sqrt(g^-1 f)), with the actual S6.20
beta2=beta3=0. The canonical clock phi and free canonical chi remain
coupled to the physical g. No change of the matter action is implicit.

Use the local coefficient interval |u|<=1/10 and 2<c<=4, with M,tau>0.
The algebraic solution of f's own zero-derivative equation is

    f0 = r(phi)^2 g,  r^3 = -beta1/beta4,  r>0.

Here r is neither the fixed action parameter c nor the rolling spatial
ratio y. The literal positive beta1 and negative beta4 make this root
and its full symmetric potential Hessian regular for every admitted
finite parameter. There is no uniform coefficient bound as c tends
to2, and c=2 is not an action in the domain.

The expansion assigns one formal derivative to each spacetime
derivative and keeps coefficient functions of phi as coefficients.
It recursively solves f's own equation about f0. It does not solve
the physical g equation for f and substitute that other equation into
the action. The exact nonlocal elimination additionally requires a
specified inverse operator and homogeneous data. No such prescription
or small remainder is established by this formal calculation alone.

## Retained action and tensor test

Modulo displayed boundaries, its zero- and two-derivative density is

    -M^2(1+r^2)R_B/2 + (1/2-3M^2 r_phi^2)X + Y/2
    -2(beta0+3r beta1).

Write B_mn=G_mn[r^2 g], with its covariant components contracted using
the physical g. The four-derivative correction is

    kappa [B_mn B^mn - (tr_g B)^2/3],
    kappa = M^4 r^3/(4 beta1).

The full traceful potential Hessian, conformal curvature identities and
constant-r physical-source Schur kernel are checked independently.
After covariant integration by parts, retain the separate operator
kappa(R_mn R^mn-R_B^2/3); do not discard it as a boundary. The remaining
quadratic scalar-tensor part satisfies A1=2F2_X and A3=0, with explicit
F and K coefficients and all derivatives of r and kappa retained.

A literal physical tensor variation proves that the retained action
has fourth-order tensor coefficient kappa, while every action in the
original quadratic scalar-tensor basis has a second-order tensor
equation on this test. At u=0, kappa=M^2 tau^2 c(c-2)/16>0. Therefore
the retained action is not exactly CD, or any quadratic-DHOST action,
in the same metric modulo boundaries and scalar-only clock changes.
This is not a ghost verdict on the untruncated parent.

## Original clock, normalization and necessary remainder size

The prescribed CD clock is theta=tau u. Its relation to the canonical
parent clock is fixed by the background:

    phi = M integral_0^(theta/tau) sqrt(kbar(v,c)) dv,
    X_phi = J(theta)^2 X_theta,  J=dphi/dtheta>0.

For a favorable leading tensor normalization, grant M_*^2=5M^2 at the
center. This does not also retune the fixed free-chi trajectory.
At theta=0, on the entire open tube 9/10<X_theta<11/10, the retained
scalar-tensor normal form has

    F2/M_*^2=-1/2,  F2_X/M_*^2=A1/M_*^2=A3/M_*^2=0.

The frozen CD target instead has

    F2/M_*^2=-X_theta/2,  F2_X/M_*^2=-1/2,
    A1=0,  A3/M_*^2=1/X_theta.

The combination Xi=A1-2F2_X transforms as Xi_theta=J^2 Xi_phi.
Consequently I=X Xi/G_T is invariant under scalar-only clock changes
and constant action normalization. In this specified scalar-tensor
normal form it is0 for the retained result and1 for CD at the center.
The independent tensor-symbol obstruction applies to the full retained
action without this normal-form qualification.
Here G_T is the coefficient of the displayed scalar-tensor normal form,
not the entire kinetic response after folding in derivatives of kappa.
Those derivatives are retained in the full fourth-order tensor equation.

If an omitted scalar-tensor coefficient remainder Delta exactly repairs
the target in this same form, it must satisfy at the center

    Delta F2_X=-M_*^2/2,  Delta A1=0,
    Delta Xi=M_*^2.

More generally, if each final normalized A1 and F2_X coefficient is
within epsilon of CD, then

    |Delta F2_X|/M_*^2 >= max(0,1/2-epsilon),
    (|Delta A1|+2|Delta F2_X|)/M_*^2 >= max(0,1-3epsilon).

These are necessary coefficient-norm bounds, not a measured actual
remainder and not a bound in every possible physical-response norm.
The target's inverse-X coefficient also means that total-derivative
grading alone is not an expansion on its X near1 tube. Higher formal
orders, resummations or a different controlled description are not
excluded; they have to supply and bound the missing structure.

## Actual background derivative defect

Let f2 be the additive order-two covariant metric correction, not the
relative perturbation h=f2/r^2. Retain the physical T component basis and
put R00=f_actual,TT-f0,TT-f2,TT. At u=0 and delta=c-2,

    f2_TT=4delta,  R00=delta^2,
    (f2_TT)_uu ->64,  R00_uu ->0,  R00_uuuu ->10752.

The fourth-derivative defect is strictly greater than10752 for every
2<c<=4; its cleared difference is an explicit polynomial with positive
coefficients in delta. For0<delta<=1/100, the first correction's second
derivative is greater than60. With varphi=phi/M, differentiating these
same physical components as functions of the clock label gives

    (f2_TT)_varphi2 ->6400/2399,
    R00_varphi4 >107520000/5755201.

The nonlinear clock conversion includes the kbar_uu term. These are not
the metric components after a spacetime-coordinate transformation to
varphi; the latter would have an additional Jacobian factor. Physical
T derivatives carry the corresponding powers of1/tau.

Thus delta alone gives neither a C2-small first metric correction nor
a C4-small retained metric defect on a fixed neighborhood in these
specified charts. On fixed compact inner intervals u=sqrt(delta)x,
the value defect instead has the consistent expansion
delta^2(1+64x^2+448x^4)+O(delta^3). Neither fact is promoted to a
norm-independent exclusion of finite-band or resummed descriptions.

## Matter transformations and remaining obligations

An explicit inverse-metric redefinition can cancel the curvature-square
term against the leading Einstein term, but it induces

    (kappa/Q)[R_mn chi^m chi^n-R_B Y/3],  Q=M^2(1+r^2),

as well as clock/Q-derivative terms. Leading-equation rewriting produces
a pure-chi contact2 kappa Y^2/(3Q^2), together with mixed clock contacts.
The full source and matter dictionary must be retained; this is not
the unchanged free M1 action in a new supposedly physical metric.

No rolling spectral gap, convergent derivative expansion, omitted
operator/loop bound, vacuum positivity, finite-gravity dispersion
remainder, global healthy completion, whole-row exclusion or original
P8 completion is claimed. A small background lapse residual is not
substituted for the open-tube off-shell coefficient test. Original P8
remains open, and no new user authority is required to continue.
