# Whole source and complete spatial ADM convention

Use P8 signature +--- and the S174 algebraic chart

    g_phys=C g_hat+(1-C)du^2/X, C=R^-1/2, X=N^-2.

The lapse and contravariant shift are unchanged; h_phys=C gamma.
The source-pinned chart has explicit inverse and covariant metric Jacobian
C^9. The complete generic-R time and spatial integrations by parts remove
lapse velocities and light-action spatial lapse gradients. This fact
belongs to the full action, not its clock polynomial. No regular unitary
chart is asserted at a zero clock gradient.

Write V=sqrt(det gamma), M=R^1/4, U=R^-3/4, Cchi=R^-1/4,
C3=R^3/4/2, and

    B=-U R_u/(2N)-I,
    Fhat=U[F+9 R_u^2/(16R N^2)]-I_u/N,
    I_N=3 U R_u R_N/(4RN), I(u,1)=0.

The primitive includes its full time and spatial boundary, not only a
coefficient evaluated on the clock. The full current bindings in the
report retain R,F, their actual fixed-profile bindings and every vacuum
constant; the normalized heavy source is exactly

    j(u,N)=physical_source(u,N^-2)/sqrt(kappa).

It is generally nonzero, and its lapse derivatives are retained. The
heavy field h and its two canonical terms are U(-n h^2/2+jh) in the
Lagrangian, with the unchanged n=mass_squared. Both matter gradients use
Cchi rather than U because the inverse physical spatial metric also
transforms. All normalized formulas factor out the same fixed kappa;
restoring it only restores the fixed canonical normalization, not a
coefficient-dependent Jacobian cancellation.

Let T=(W0-Ni Wi)/N and K=tr Khat. The original source is a one-form along
du, with normal component

    Snormal=(R-1)(K-3Hclock(u)/N), Hclock=4u/(1+u^2).

Hclock is a fixed action coefficient, not the varied Hhat. Set r=R-1 and
c=-3Hclock r/N. The trace/temporal Lagrangian per N V is

    Ltrace=-(M/3)K^2+B K+Fhat+U(T-rK-c)^2/2.

Keep the five traceless metric velocities, all three electric velocities,
the magnetic and spatial vector mass terms, both matter velocities and
gradients, and C3 times the intrinsic spatial scalar curvature. Their
independent Legendre calculations are part of the packet.

The complete electric Legendre transform contains pi_W^i partial_i W0.
Its spatial integration by parts gives -N T partial_i pi_W^i and the
retained shift terms. Thus it removes lapse derivatives from the
auxiliary Hamiltonian without setting the Gauss density to zero. The
source square is neither aligned nor removed before variation.

