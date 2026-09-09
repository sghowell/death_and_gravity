# Fresh global canonical chart and explicit domain

All coefficients here are rederived from S6.98's full original phase,
not from either off-clock solution. On |u|<=1/4 the gamma pivot
Lambda is separated from zero. First make the exact canonical
exchange

    b=-pi_v/(2*a^3*q), P_b=2*a^3*q*v,
    x=(b,chi), p=(P_b,pi_chi).

It is time dependent at fixed comoving momentum. The derivative
of its matrix is retained when computing E'+EM, so the moving
boundary term is not lost.

The leading velocity-coordinate coefficient in the full gamma
Lagrangian is the symmetric B_q=diag(-2*Theta/Lambda,0).
Shift p_tilde=p-a^3*q*B_q*x, and use the packet variables

    Y=(k*a^(3/2)*x, a^(-3/2)*p_tilde), k=|xi|>0.

The explicit map E from density phase to Y satisfies
E Omega E^T=k Omega. Accordingly E/sqrt(k), NOT E, is symplectic.
Its exact generator is computed as (E'+EM)E^-1, including q'=-2Hq.
It has the complete finite Laurent form

    Y'=(k*J0+L0+L1/k)Y.

The possible inverse-square and inverse-cube terms vanish identically
for this GLOBAL clock; their vanishing is checked, not assumed.
Every retained coefficient is Hamiltonian relative to constant Omega.

Put V=[1,0;w/Lambda,1], kappa_c=2*(J+delta_J)/Lambda^2,
D=diag(kappa_c,1), omega=diag(c/a,1/a), c^2=J/(J+delta_J).
The leading Hamiltonian blocks are

    H_pp=K^-1, K=V^-T D V^-1,
    H_qq=V^-T D omega^2 V^-1, H_pq=H_qp=0.

The complete moving gradient identity is substituted from the actual
global background to verify these formulas. Both leading blocks
are positive. B0=V^-T D omega V^-1 is positive and
R0=-i B0 solves the leading Riccati equation. Configuration evolution
on this graph has frequencies -i*k*c/a and -i*k/a.

At u=0, kappa_c=243/20, c=sqrt(1199/1215), a=1 and
w/Lambda=-1/10. Thus

    B0=[(243/20)*c+1/100, 1/10; 1/10,1].

The exact order-zero Hamiltonian and first derivative R0' vanish
at this time. R1(0)=0, but R1'(0) and R2(0) DO NOT vanish.
Their explicit native values and the full recursion residuals are
included in the report. Dropping the time derivative gives a
different second correction.

For global microlocal propagation we also need separation on every
finite strip. Write J*h^2=P(u)/[800*(1+u^2)^12]. The literal numerator
P has positive even coefficients, checked entrywise. If |u|<=R,
R>=1/4, then J*h^2<=P(R)/800=:B_R. Let A_R=(1+R^2)^2.
Consequently

    1-c^2 >= 4*epsilon/(B_R+4*epsilon)=:d_R>0,
    (1-c)/a >= d_R/(2*A_R),
    2*c/a > 493/(250*A_R).

These explicit finite positive bounds cover same-sign and opposite-
sign frequency gaps. They do not assert an asymptotic-infinite-time
gap. The existing original-phase evolution is used below the
large-frequency chart threshold, including every gamma velocity pole.
