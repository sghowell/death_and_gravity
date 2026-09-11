# Vary the covariant action before fixing the lapse

Use the P8 physical (+---) metric ds^2=N^2dt^2-a^2dx^2, unit
clock u=t, and v=partial_t log a. The dimensionless kinetic
invariant is X=N^-2. Direct metric connection and contractions give

    R=-6[vdot/N^2-v Ndot/N^3+2v^2/N^2],
    box u=-Ndot/N^3+3v/N^2,
    Z=u^mu u_;mu nu u^nu=-Ndot/N^5,
    L3=Ndot^2/N^8-3v Ndot/N^7,
    L4=Ndot^2/N^8, L5=Ndot^2/N^10.

K=A1=A2 vanish identically in the frozen target. For its first
metric variation use the exact clock jets

    F=F0(t)+FX(t)(X-1)+O((X-1)^2),
    F2=-1/2-(X-1)/(2h(t))+O((X-1)^2),
    A_i=A_i(t,1)+O(X-1).

The O terms do not affect the first background variation: L3 is
already first order in lapse velocity and L4,L5 are second order.
These are local jet representatives, not replacements of the full
analytic target away from the clock.

For the density L=N a^3 times the covariant kernel, first compute

    rho_Euler=-(partial_N L-d_t partial_Ndot L)/a^3,
    P_Euler=(partial_b L-d_t partial_v L+d_t^2 partial_vdot L)/(3a^3),
    b=log a,

and only then put N=1,Ndot=Nddot=0. The Einstein part gives
rho_E=-3H^2, P_E=2Hdot+3H^2. The scalar gives rho_F=2FX-F0,
P_F=F0. Relative to that fixed Einstein part,

    rho_curvature=6(Hdot+2H^2)/h, P_curvature=0,
    rho_A3=-3[A3 Hdot+Adot3 H+3A3 H^2], P_A3=0.

A4,A5 give zero first background Euler source. That does not mean
they can be deleted from the action or perturbation constraints.
All quantities here are normalized by kappa, the overall physical
action multiplier; they are Euler-source components of one action,
not independently conserved physical matter tensors.

The zero clock values of the extra curvature and A3 kernels do not
imply zero metric variations. Setting X=1,Ndot=0 in the action
before variation would discard precisely the displayed terms.
No division by H or a crossing constraint determinant is used.
