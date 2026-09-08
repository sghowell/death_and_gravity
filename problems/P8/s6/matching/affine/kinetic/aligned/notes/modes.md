# Actual canonical vector frequencies on the rolling background

Let a=(1+u²)^2, H=a'/a, q=k_com²/a², zeta>0, and r=zeta*q.
Primes hold comoving momentum fixed, so q'=-2Hq. These are the actual
decoupled quadratic equations obtained after the full constraints,
not an isolated flat-space mass estimate.

For either transverse coordinate component T_i,

    L_T=a*[zeta*T_i'²-(zeta*q+1)*T_i²]/2,
    g_T²=a*zeta, v_T=g_T*T_i.

For the nonzero-momentum longitudinal potential,

    L_L=a³*[zeta*q*sigma'²/(1+zeta*q)-q*sigma²]/2,
    g_L²=a³*zeta*q/(1+zeta*q), v_L=g_L*sigma.

Both canonical Euler equations are v''+Omega²*v=0 with

    Omega²=q+1/zeta-U, U=g''/g,
    U_T=H'/2+H²/4,
    U_L=alpha(r)*H'+beta(r)*H²,
    alpha=(1+3r)/(2(1+r)),
    beta=(1-2r+9r²)/(4(1+r)²).

In particular g_L'/g_L=H*alpha and
beta=alpha²-2r*alpha_r. Holding physical q fixed would change this
answer; the tests detect that error and omission of the normalization.

## Uniform all-time bound

For r>=0, 1/2<=alpha<3/2 and 0<beta<9/4. The latter follows from

    1-2r+9r²=9(r-1/9)²+8/9,
    9/4-beta=(2+5r)/(1+r)².

Exact polynomial identities give H²<=4 and |H'|<=4:

    4-H²=4(u²-1)²/(1+u²)²,
    4-H'=4u²(u²+3)/(1+u²)²,
    4+H'=4(u^4+u²+2)/(1+u²)².

Thus |U_T|<=3 and |U_L|<=15, yielding

    Omega²>=q+1/zeta-15>=q+1985, 0<zeta<=1/2000,

at every finite real u and every positive comoving momentum. At zero
momentum g_L=0 is not an invertible normalization; instead there are
three homogeneous coordinate-vector components governed by the same
transverse weight and U_T. This separate chart obeys the same bound.
The r->0 correction agrees with U_T, which is a consistency control,
not permission to invert the singular longitudinal map at k=0.

With physical units restored, zeta=zeta_physical/(M²*tau²) and

    Omega_physical²>=k_physical²+M²/zeta_physical-15/tau².

For example M²=3, tau=2, zeta_physical=1/1000 gives zeta=1/12000.
At u=1/2 and k_com²=4375/256, q=7 and the physical lower bound is
2998. Unit normalization is checked explicitly, not inferred by M=1.

This frequency floor is a property of a time-dependent quadratic
initial-value equation. It is not a stationary spectral gap, a heavy
vacuum S-matrix threshold, a nonlinear degree-count theorem or an EFT
cutoff. A retarded inverse and its source/state-dependent errors must
be established separately before any heavy-elimination claim.
