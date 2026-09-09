# Literal primitive-free rational evolution

Let the homogeneous Hamiltonian be

    Hcal=a(p-b)^2-Nf0+m ell^2,

with the actual original boundary primitive I. Define d=b+I,
V=Nf0+I_phi and y=p-b. Neither d nor V contains a primitive.
Use the literal endpoint I_N=-I_s(u,1/N)/N^2; its time derivative
is retained. At fixed original p,ell,

    F=Hcal_N=a_N y^2-2a y b_N-V_N+(I_N)_u+m_N ell^2,
    h=F_N|y-b_N F_y.

From the original canonical background equations,

    y'=-a y^2+V-d_u+m ell^2-b_N N',
    ell'=-2a ell y.

I_phi cancels between p' and b'; it is not set to zero in either.
This cancellation is valid for the same fixed-basepoint primitive,
with commuting partial derivatives on its analytic branch.

## Rational variables

Set z=e y, w=e ell and

    hbg=(1+u^2)^3, D=1+N^2(hbg-1), e^4=N^2 hbg/D,
    Ln=e_N/e=1/(2ND),
    Lt=e_u/e=3u(1-N^2)/(2(1+u^2)D).

Write a=e a0, m=e m0, d=d0/e, V=V0/e and I_N=I0/e.
The code derives d0,V0,I0 from the literal coefficients, not from a
Taylor replacement. In particular,

    a0=-3N/4, m0=D/(2N hbg),
    d0=-3u(N^2-1)/(N(1+u^2)D),
    I0=-9u(N^2-1)/(N^2(1+u^2)D^2).

V0 is an exact rational polynomial quotient retained in the report.
Define

    b0=(d0)_N-Ln d0-I0,
    A=(a0)_N+Ln a0,
    M=(m0)_N+Ln m0=[2N^2(hbg-1)-1]/(4N^2 hbg),
    B=-2a0 b0,
    C=-(V0)_N+Ln V0+(I0)_u-Lt I0,
    f=eF=A z^2+Bz+M w^2+C,
    Kz=Ln z-b0, Kw=Ln w,
    Z0=-a0 z^2+V0-(d0)_u+Lt d0+m0 w^2+Lt z,
    W0=(Lt-2a0 z)w,
    PIV=f_N+f_z Kz+f_w Kw,
    FORCE=f_u+f_z Z0+f_w W0.

The actual lapse Hessian is (PIV-Ln f)/e. Only on f=0 may
one write h=PIV/e. Keeping that off-constraint term prevents
confusion between fixed-phase derivatives and differentiation
along the family. The evolution is

    N'=-FORCE/PIV, z'=Z0+Kz N', w'=W0+Kw N',
    Hhat=-Nz/2.

The identities follow from the original equations, not an
independently prescribed background. The solved lapse velocity
makes the actual derivative of f vanish.

## Actual characteristic

Since alpha=f_z/3 and r_N=-eM, put Bbar=-f_z/(3M), so the
principal B of reduction.md is Bbar/e. Define

    G0=(Bbar)_u+(Bbar)_z Z0-(Lt+Hhat)Bbar-m0 Bbar^2-Nw^2,
    GN=(Bbar)_N+(Bbar)_z Kz-Ln Bbar.

Bbar is independent of w. Then

    e G_clock=G0-GN FORCE/PIV,
    c_clock^2=-4 hbg M^2 (G0 PIV-GN FORCE)/(D PIV^2).

All expressions are rational. Only after the fixed-phase lapse
and time derivatives above are formed, impose

    w^2=-(A z^2+Bz+C)/M.

The resulting rational function is exact on the whole constraint
branch, not a finite Taylor jet. Its central value reproduces S6.83
as a rational function of independent N. Its numerator and
denominator have 1457 and 1387 nonzero monomials; the report pins
their exact sorted monomial/rational-coefficient hashes as well as
all intermediate rational-system hashes. Rebuilding performs the
complete arithmetic, not a lookup of those hashes.
