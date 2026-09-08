# Covariant quadratic pole and independent geometric derivation

Use epsilon=(4-d)/2 and Delta=-nabla^2. In units of
1/(32 pi^2 epsilon), the zero-derivative quadratic density is

    Q0=-Y_mn Y^mn/8-(tr Y)^2/16.

The nine two-derivative invariants, including curvature times Y^2,
are written literally in invariants.SECOND and divided by 48m^2.
Ordered derivatives are not commuted before evaluation. An independent
first-derivative basis from Garcia-Recio--Salcedo (114), with its
epsilon sign converted, is in kernel.independent_second. Both
independent temporal and spatial mass-profile Euler variations
of the difference vanish, with arbitrary H and physical momentum.

## Why the four-derivative term is a scalar determinant

Write A=1+Y/m^2 and use the Stueckelberg shift W=B+dvarphi/m.
Lorenz gauge gives a vector block Delta_H+m^2 A, a scalar block
G=-nabla_mu A^mu_nu nabla^nu, an off-diagonal block m A nabla,
and a mass-independent scalar ghost. The block determinant
expansion has only even off-diagonal insertions. Its second
term carries m^2 and has at most two background derivatives in
its UV residue; its fourth carries m^4 and has none. Higher
terms are UV finite. The vector and ghost four-derivative
residues have no A dependence. Consequently every four-derivative
Y-dependent residue comes from (1/2)Tr log G.

Set

    g_tilde=sqrt(det A) A^(-1) g.

In four dimensions sqrt(g_tilde) g_tilde^(-1)=sqrt(g) A g^(-1).
The scalar quadratic action is the ordinary massless scalar action
for g_tilde. The ultralocal field-measure change has a scaleless
Jacobian in dimensional regularization. The universal minimal
scalar coefficient gives

    Q4=-[lambda^2] {sqrt(g_tilde(lambda))/sqrt(g) a4sc(g_tilde(lambda))},
    A(lambda)=1+lambda Y/m^2,
    a4sc=R^2/72-Ricci^2/180+Riemann^2/180.

Only the covariant total Laplacian is omitted. This definition
retains the full curvature basis; it does not erase a
four-dimensional Gauss-Bonnet term before a future dimensional
continuation. It is a local quadratic functional of arbitrary Y.
The executable component reduction is the FLRW isotropic mass
sector needed by the actual first lapse insertion.

## Two independent geometric representations

In the Euclidean original frame, e0=partial_u and ei=a^(-1)partial_i,

    nabla_e0 ei=0,
    nabla_ei e0=H ei,
    nabla_ei ej=-H delta_ij e0.

Ordered derivatives of Y include connection terms for every
inner derivative index as well as both tensor indices. Physical
momentum obeys k'=-Hk. The Riemann tensor has sectional entries
R_0i0i=-(H'+H^2), R_ijij=-H^2. Ricci contraction, trace
Laplacian and covariant commutator identities check this convention.

Independently, geometry.py constructs the coordinate Christoffel
symbols and full Riemann tensor of

    ds_tilde^2=N_tilde^2 du^2+A_tilde^2 dvec(x)^2,

allowing dependence on u and one spatial coordinate. With
t=Yt/m^2,s=Ys/m^2,

    N_tilde^2=(1+s)^(3/2)(1+t)^(-1/2),
    A_tilde^2=a^2 (1+t)^(1/2)(1+s)^(1/2).

The full scalar heat density is expanded to second order in t,s
and their coordinate jets. Only afterwards is its coefficient
polarized into opposite spatial Fourier legs. The original a
can be set to one at the evaluation point if every coordinate
spatial jet is rescaled by the same a; its time derivative H
is retained. This is not a homogeneous truncation.

## Independent checks and boundary conditions

For H=0 and time jets (+/- i omega)^j, all three coefficients
equal minus one half of the frozen S6.48 Feynman parameter
integral. The factor one half is fixed by the prefactors
1/(32 pi^2 epsilon) here and -1/(64 pi^2 epsilon) there.

For Y=m^2 f g, g_tilde=(1+f)g. Four-dimensional Maxwell conformal
invariance maps the full operator to ordinary constant-mass Proca
on that conformal metric. Modulo compact-variation boundaries,

    Q0=-3m^4 f^2/2,
    Q2=-3m^2 (grad f)^2/4,
    Q4=-[9(Box f)^2+6R f Box f+3R(grad f)^2]/72.

Here Box f=f''+3Hf'-k^2 f and (grad f)^2=f'^2+k^2f^2 in
the opposite-momentum bilinear action. Weighted Euler variations
use the formal adjoint -D-3H and vanish exactly for all three
differences. For constant independent t,s, the scalar heat
density instead scales by [(1+t)/(1+s)]^(3/2); its quadratic
coefficient checks the nonconformal curvature dependence.

These are local UV residues. Neither their signs nor their
polynomial high-frequency behavior is a finite causal feedback,
stability, exact higher-derivative-mode or cutoff verdict.
