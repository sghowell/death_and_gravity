# Complete finite-width costs and exact null margin

Work in dimensionless affine x=lambda/tau. For either clamped cubic
piece of length l, the squared derivative norms of orders zero,
one and two are

    (n0,n1,n2)=(13l/35,6/(5l),12/l^3).                   (12)

Both pieces are integrated directly. Use l=1/1000 for the past and
l=2 for the future. Their spatial profile has squared first derivative
norm 3 in the dimensionless spatial coordinate z/tau.

On a segment with lo<=A<=hi and derivative bounds d1,...,d4,
Cauchy-Schwarz and (4) give the dimensionless majorants

    BP=3[hi^2 n2+4d1^2 n1+4(d1^2/lo+d2)^2 n0],
    BM=2[hi^2 n1+4d1^2 n0]/lo^4,
    CQ=(5/9)BP+BM,
    CW=2n1+8d1^2 n0/lo^2+3n0/lo^4.                     (13)

The inequalities (u+v+w)^2<=3(u^2+v^2+w^2) and
(u+v)^2<=2(u^2+v^2) hold pointwise; no cancellation or sampled-time
positivity is required. CQ retains the mixed spatial term and CW
retains its spatial state term. Bounding both signs of the actual
reference (7) gives

    B0=4d1^2 d2/lo,
    Bbeta=48d1^2 d2/lo+84d2^2+72d1d3+12hi d4.            (14)

The reference costs per segment are n0 B0/360 and n0 Bbeta/360,
because kappa hbar/(2880pi^2)=delta tau^2/360. No photon coefficient,
trace-free simplification or beta_S=0 convention is substituted.

Insert the past and future tuples in (P)-(F). Sum CQ plus the zero
scheme reference to obtain C0, sum finite-scheme references to get
Cbeta, and sum CW and n0 for Cstate and Csource respectively.
All quantities in (12)-(14), both segment contributions and their
totals are reconstructed in native SymPy and a separate Fraction-only
implementation. The exact totals are displayed in the formulation.

In particular C0>Cbeta>0. Thus

    delta(C0+|beta_S|Cbeta)
       <=delta(1+|beta_S|)C0<=10^-12 C0.

The complete maximum index upper bound is

    6/5+104/7875+10^-12 C0+10^-6 Cstate+Csource/10
      =6843647823531587849/4784062500000000000,

strictly below 3/2 and hence below the outgoing threshold 2. This
leaves the exact margin in (2), greater than one half. All terms
are nonnegative, so the same margin holds throughout the full
rectangular budget domain. No optimum is claimed for the cubic
sampler, derivative majorants, future caps or the final constants.

Macroscopic scaling means delta is proportional to kappa hbar/tau^2
and the remaining budgets are explicit dimensionless field/source
conditions. For fixed finite beta_S, the quantum requirement imposes
only a finite Planck-to-affine-time separation. It is not a measurement
of the state cap or finite renormalization in an observed cosmology.
