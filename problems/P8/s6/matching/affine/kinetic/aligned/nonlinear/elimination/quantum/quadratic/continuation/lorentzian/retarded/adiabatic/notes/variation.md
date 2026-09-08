# Vary the physical mode readout before dimensional integration

Use the same mass source a_m=1+alpha n,b_m=1+beta n at a fixed metric.
In D spatial dimensions, for physical q and omega^2=q+m^2,

    g_T^2=a^(D-2), Omega_T^2=q+m^2 b_m,
    g_L^2=a^D m^2 a_m q/(q+m^2 a_m),
    Omega_L^2=b_m(q/a_m+m^2).

At zero source let z=q/omega^2 and lambda=omega'/omega=-H z.
The background rates are d_T=(D-2)H/2 and d_L=d_T+H z.
Direct differentiation of the physical Hamiltonian gives

    r_T=delta log Omega_T=beta(1-z)n/2,
    r_L=delta log Omega_L=(beta-alpha z)n/2,
    r_g,T=0, r_g,L=delta log g_L=alpha z n/2,
    delta d=r_g', delta lambda=r'.

These derivatives commute with time evolution at the fixed metric.
The time derivative acts on q through q'=-2Hq, or equivalently
z'=-2H z(1-z), and on all source, alpha, beta and Hubble jets.

For the exact instantaneous mass insertion, in the moving canonical
variables v=g w and p=Pi/g, J=(A|p|^2+B Omega^2|v|^2)/2.
At nonzero source the weights are

    A_T=0, B_T=-m^2 beta/Omega_T^2,
    A_L=alpha q/[a_m(q+m^2 a_m)], B_L=-beta/b_m.

Their variations at zero source are delta B_T=-2r_T B_T,
delta A_L=-alpha^2 z(2-z)n and delta B_L=beta^2 n.
These are derived independently, not obtained by dropping the
S6.64 temporal-constraint contact. The moving mode/readout
derivatives together represent that full physical variation.

## Local reference expansion

Set U=d'+d^2 and W=omega(1+P2/omega^2+P4/omega^4), with

    P2=-U/2-lambda'/4+lambda^2/8,
    P4=-P2^2/2-P2''/4+5lambda P2'/4
       +(lambda'/2-3lambda^2/2)P2.

The linearized recurrence retains delta U=(delta d)'+2d delta d,
delta P2 and delta P4, including all alpha/beta/source derivatives.
The first two coefficients of the varied exact Riccati residual
vanish identically in arbitrary D; this independently checks
the differentiated recurrence.

For c1=d+lambda/2 and b2=P2'-2lambda P2, the three local readout
coefficients multiplying omega^(1-2j)/4 are

    C0=A+B,
    C1=(A-B)P2+A c1^2,
    C2=(A-B)P4+B P2^2+A(c1 b2-c1^2 P2).

The actual varied coefficient is delta Cj+(1-2j)r Cj.
The transverse multiplicity remains D-1 until after integration.
Only the local expansion is being varied here. The exact selected
mode response minus this expansion still needs its uniform
integrability and quantitative response estimates.
