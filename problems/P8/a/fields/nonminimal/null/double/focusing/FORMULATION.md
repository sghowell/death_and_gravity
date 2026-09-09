# P8-A.23 — null incompleteness from the physical scalar double-null QEI

This is a conditional null-geodesic theorem for a stipulated GLOBAL
spatially flat FLRW spacetime. Its energy input is the actual physical
conformal scalar double-null QEI of A.22, with a one-sided relative
Wick-square upper cap. It does not assume a homogeneous quantum state,
a single-null-line QEI, pointwise NEC or pointwise SEC. It does not
assign original P8 closure.

## Geometry, field and affine normalization

Let M=I_t times R^3 with g=dt^2-a(t)^2 dX^2, smooth a>0, the A.20
FK curvature convention, and the actual semiclassical equation

    G_FK+Lambda g_FK=-kappa(T_scalar+T_other), kappa>0.    (1)

The scalar is one real massless field with xi=1/6 in any Hadamard
target, in the fixed scalar stress scheme of A.20 with finite real
beta_S. Lambda is arbitrary and cancels from the null contraction.
Normalize a(0)=1, choose tau>0, and assume

    H(0)<=-2/tau.

Use the future affine clock lambda(t)=integral_0^t a(s)ds, so the
fixed null vector is K=a^-1 partial_t+a^-2 partial_z. Its initial
proper energy is one. Set x=lambda/tau and A(x)=a(t(lambda)).
Derivatives in the following two caps are with respect to x,
NOT proper time. There is an actual past interval [-1/1000,0],
with smooth endpoint neighborhoods, on which

    9/10<=A<=2,
    |(A',A'',A''',A'''')|<=(3,16,1000,2000000).           (P)

IF the future affine endpoint exceeds 2tau, require throughout
x in [0,2] the separate conditional caps

    1/2<=A<=2,
    |(A',A'',A''',A'''')|<=(50,10000,2000000,400000000).  (F)

Neither (F) nor a future state cap is inferred from (P).

On the fixed coordinate plane (t,0,0,z) with |z|<=tau and the above
past affine interval, assume w=<:Phi^2:conf><=Phi_*^2, Phi_*^2>=0.
The same cap is required on the future plane sample through lambda=2tau
IF that extension exists. For an additional source, require
T_other,KK>=ell on the same finite plane sample, where ell is constant.
For the pure scalar theory this extra source is absent and ell=0.
The cap is not a bound on |w|, an operator norm or a momentum cutoff.

Define

    delta=kappa hbar/(8pi^2 tau^2),
    zeta=kappa Phi_*^2/3,
    sigma=tau^2 kappa max(-ell,0).

## Theorem

Under the displayed hypotheses and

    delta(1+|beta_S|)<=10^-12,
    zeta<=10^-6,
    sigma<=1/10,                                         (Q)

the future null affine endpoint is at most 2tau. In particular the
stipulated global spacetime is future null geodesically incomplete.
All initially unit-proper-energy null geodesics from t=0 have the
same affine duration by FLRW symmetry. Reversing the contracting
orientation gives the corresponding past-directed statement.

The proof uses a round initial sphere of physical radius tau. Its
outgoing null expansion is theta0=2(H0+1/tau)<=-2/tau. With a clamped
future affine sampler on [0,2tau], the exact sourced index bound is

    tau I <= 6/5+104/7875
             +delta(C0+|beta_S| Cbeta)+zeta Cstate+sigma Csource,

    C0=110874110523310814/1366875,
    Cbeta=3737503266601/65625,
    Cstate=23677459436/382725,
    Csource=26013/35000.

Its worst-case margin satisfies

    2-tau I >=2724477176468412151/4784062500000000000
              >1/2.                                    (2)

The outgoing Jacobi scale remains positive on any actually existing
finite affine segment. Its exact square identity instead requires
I+theta0>=0. Equation (2) is incompatible with that identity if an
extension beyond 2tau exists. This proves the stated endpoint bound.

## Physical and logical boundaries

The plane-to-curvature step uses (1): TOTAL stress is geometric.
It does not require the scalar's state or stress to be homogeneous,
and does not shrink the second sampling width to zero. The finite
width tau and all its derivative costs are retained. The full scalar
reference, including beta_S, is retained. For an additional source,
its null lower bound is an explicit additional energy premise.

A smooth future-complete geometry satisfies all the geometric caps
and the initial contraction. Its conformal vacuum also satisfies
w=0, but the pair is not asserted to solve the allowed small-source
SEE. This control shows that the geometric assumptions alone do
not already force the conclusion. Separately, the actual thermal
scalar SEE branch of A.20 at delta=10^-14 realizes (P), the initial
data and the stronger past field budget. Its earlier affine endpoint
does not establish (F) or the conditional future state cap, and is
not a new QEI-only endpoint mechanism.

The budgets are explicit macroscopic-scale inequalities for fixed
finite beta_S, not optimized constants or an observed-universe fit.
No conclusion is made about every larger extension, curvature
blowup, fundamental EFT endpoint control or general non-FLRW geometry.
The exact algebra and written analytic proof are source-pinned, not
proof-assistant formalized or independently peer reviewed.
