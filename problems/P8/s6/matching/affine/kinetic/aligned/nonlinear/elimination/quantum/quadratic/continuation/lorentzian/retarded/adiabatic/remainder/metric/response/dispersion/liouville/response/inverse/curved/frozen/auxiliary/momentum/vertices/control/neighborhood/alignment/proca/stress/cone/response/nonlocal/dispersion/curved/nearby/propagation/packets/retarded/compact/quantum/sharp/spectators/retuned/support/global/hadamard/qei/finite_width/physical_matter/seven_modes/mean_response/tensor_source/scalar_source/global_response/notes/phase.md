# An actual global scalar phase estimate

Put t=1+u^2, a=t^2. Let Z=(v,chi,p_v,p_chi) be the
S6.103 canonical density phase. S6.106 derives the exact
symplectic bridge Z_natural=B(u)Z, including BOTH momentum
shifts and the generating time boundary. The actual added
anchor covariance is eta_S b I4 in that natural phase.

For estimating the same free evolution use NONCANONICAL
coordinates W=T Z with

    T=diag(1,1,(2/5)t^-3,(2/5)t^-3),
    K=T M_old T^-1+T' T^-1,
    R=K+diag(0,0,6u/t,6u/t).

For u>=0 the removed diagonal is nonpositive. The upper
Dini derivative of the infinity norm is at most the
maximum absolute row sum of R times the norm; replacing
that maximum by the sum of all absolute entries is valid.
No positive damping or off-diagonal coefficient is dropped.

Each rational entry is bounded on 1<=k<=2 by

    C |u|^delta/t^r, delta in {0,1}.

The denominator has nonnegative even-time coefficients
and positive constant. After replacing numerator coefficients
by their absolute values and k powers by their value at 2,
the code checks that C times that denominator minus
t^r times the numerator divided by u^delta has nonnegative
coefficients. Unsupported proof structures are rejected.
For odd entries r>=2; for even entries r>=1.

Odd half-line integrals are 1/[2(r-1)]. Even integrals use
the exact integer recurrence and pi<22/7. Half-integer
integrals used below have base I_(3/2)=1 and recurrence
I_s=(2s-3)I_(s-1)/(2s-2). These are inequalities or identities
on whole half-lines, not quadrature at sample points.

The sum of entry integrals is strictly below 59. Also
||T(0)B(0)^-1||_infinity=1 exactly. Thus every row of
the natural-anchor to W propagator has absolute row sum
at most Q=3^59, since exp(1)<3.

Retain the negative damping for the two momentum rows.
Each bottom row of R has absolute row sum <=C_j/t;
the exact C_j satisfy 1+C_j<40. Its integrating factor is t^3:

    |W_p(u)| <= Q t^-3[1+C_j integral_0^u t(s)^2 ds].
    integral_0^u t(s)^2 ds=u+2u^3/3+u^5/5
                              <=u t^2<=t^(5/2).

Since 1<=t^(5/2), |W_p|<=40Q/sqrt(t).
Consequently every row of the EXTRA weighted propagator
(v,chi,sqrt(t)W_pv,sqrt(t)W_pchi) is bounded by 40Q.
Its added covariance has every entry bounded by
C_cov eta_S, C_cov=(40Q)^2<10^60, after integrating b.

These weights do not redefine the state or its CCR.
The actual form is Omega_W=(2/5)Omega/t^3, and
K Omega_W+Omega_W K^T=Omega_W' exactly. The source report
checks that identity instead of imposing a constant CCR.
The selected addition's reflection, not reflection of the
reference state, extends all estimates to u<=0.
