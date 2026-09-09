# Baseline local evolution and the original equations

This note proves the coarse 10^-30 interval. The exact-jet argument
in enlarged-domain.md extends this same solution to |u|<=10^-7.

## Literal canonical reduction

Set the homogeneous source-free Proca field and its canonical momentum
to zero. Its linear homogeneous Proca equation and temporal constraint
then vanish identically; no clock-dependent vector source remains in
the S6.81 candidate. Homogeneous isotropic initial geometry has zero
spatial curvature and trace-free momentum, preserved by spatial
homogeneity and rotation invariance.

Let R be the hat scale, PR its canonical momentum and C the conserved
free-matter momentum. The reduced Hamiltonian is exactly

    Hcan=R^3 Hcal(u,N,p,ell),
    p=PR/(3R^2), ell=C/R^3,
    Hcal=N[(p-b)^2/(4a)-f0+ell^2/(2U)].

Differentiate before any constraint substitution. Hamilton's equations
give, independently checked from this literal expression,

    R'/R=Hcal_p/3,
    p'=-Hcal+ell Hcal_ell,
    ell'=-ell Hcal_p,
    chi'=N ell/U.

The lapse equation is Hcal_N=0. The nonzero lapse Hessian and the
implicit root from complex-domain.md allow its elimination.
Differentiating the reduced Hamiltonian adds no spurious root-derivative
term: it is multiplied by Hcal_N=0. Conversely the reduced solution
with that root obeys the original homogeneous canonical equations.

The physical metric is dt=N du, a=eomega R and U=eomega^3.
Consequently a^3 dchi/dt=ell R^3=C, exactly. Choosing any additive
initial chi value does not change these free-matter equations.
The isotropic trace equation is the spatial metric equation, and the
lapse constraint is the temporal metric equation. Momentum constraints
and trace-free equations vanish by the same homogeneous/isotropic
symmetry. The covariant diffeomorphism Noether identity recovers the
clock equation from these metric, matter and vector equations, because
phi'=1 is nonzero. No independent lapse acceleration is imposed.

## A uniform two-variable Picard contraction

For every central shift d in [0,10^-6], let ell0 be the positive
constraint-compatible central charge. Insert the newly proved
holomorphic lapse root into

    G(p,ell)=(-Hcal+ell Hcal_ell, -ell Hcal_p).

The root exists jointly on |u|,|p|,|ell-ell0|<=r=10^-24.
The original coefficient-domain Cauchy radii still apply at every
root value, so |Hcal_p|,|Hcal_ell|<=80000, |Hcal|<10000
and |ell|<1/4 give both components of G bounded by M=50000.
These are derivatives at fixed N before substitution; no derivative
of the implicit root is needed for this bound.

On the phase half-polydisc, each of the two phase derivatives of
each component is at most M/(r/2), by Cauchy's formula for the
composed holomorphic G. In the maximum norm its row Lipschitz
bound is therefore 4M/r. Take T=10^-30. The exact checks give

    M*T=5*10^-26<r/2,
    (4M/r)*T=1/5<1.

On holomorphic maps of the closed complex time disc |u|<=T into
the phase half-polydisc, the integral map

    (p,ell)(u)=(0,ell0)+integral_0^u G(s,p(s),ell(s)) ds

is strictly invariant and contractive. Holomorphic integration is
path-independent inside the disc; the straight radial segment gives
the same sup-norm bound. Completeness, uniform convergence and
conjugation symmetry yield a unique real analytic solution on the
real interval. Local uniqueness also holds among C1 solutions in
this domain by the same locally Lipschitz vector field.

Reconstruct log R by its exact equation and R(0)=1. Then

    |log R| <=(80000/3)*T<1/2,

so 1/2<R<2 on the real interval. The exact equations imply
(ell R^3)'=0 and keep C nonzero. The lapse is positive, remains
inside its recentered root disc, and X=N^-2 remains in the original
covariant clock tube. Thus the construction uses the original
smooth branch throughout, not a Taylor-polynomial substitute.

## Parity and a genuine physical local bounce

The actual coefficient functions a,U,f0 are even in u and b is odd.
The literal boundary primitive is odd and its time derivative is
even. Hence Hcal(-u,N,-p,ell)=Hcal(u,N,p,ell). Constraint-root
uniqueness preserves this symmetry. Evolution uniqueness implies
p is odd and ell,N,R are even; the conformal factor and physical
a are even too. In particular N'(0)=R'(0)=ell'(0)=0 and H(0)=0.

The twice-differentiated constraint computation in acceleration.md
determines the actual N''(0), including I_uuu. On the full central
family it gives dH/dt(0)>39999/10000 and a(0)=sqrt(N0)>=1.
Continuity makes a strict physical local minimum and the nearby
contracting/expanding signs follow on some neighborhood of zero.
The quantitative interval |u|<=T is a solution-existence bound:
we do not silently claim these signs or a>=a(0) on every point of
that interval without an additional variation estimate.

This realizes the S6.83 central cone data on actual local classical
solutions, including the superluminal central control. It neither
proves a global nearby bounce nor transfers a quantum response or
an interacting propagation window to the new background.
