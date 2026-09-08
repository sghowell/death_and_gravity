# Exact Liouville mode reduction

## Positive mode clocks and canonical variables

Let k>0 be comoving momentum, q=k^2/a^2, and use the physical
constrained Hamiltonians already derived in S6.67. Their kinetic
coefficients and frequencies are

    g_T^2=a/N, Omega_T^2=N^2(q+m^2 bm),
    g_L^2=a^3 m^2 am q/[N(q+m^2 am)],
    Omega_L^2=N^2 bm(q/am+m^2).

For either sector choose d sigma=f dt with f>0. A physical oscillator
coordinate Q has action g^2(Qdot^2-Omega^2 Q^2)/2. In sigma its
kinetic coefficient is C^2=g^2 f and its frequency is Omega/f.
Thus v_sigma=C Q=sqrt(f) v_t is canonical. Integration by parts
of its cross term gives

    v_sigma''+[Omega^2/f^2-C''/C]v_sigma=0.

The momentum entering the physical readout is
p_sigma=v_sigma'-(C'/C)v_sigma=p_t/sqrt(f).
The map diag(sqrt(f),1/sqrt(f)) is symplectic; the two-mode
Wronskian is unchanged. In particular,

    (A_read |p_t|^2+B_read Omega^2 |v_t|^2)/2
       =f(A_read |p_sigma|^2
             +B_read(Omega/f)^2 |v_sigma|^2)/2.

The prefactor f is part of the physical current. It must not be
discarded when varying the stress. A direct second derivative of
C=g sqrt(f) verifies the time-change terms

    C''/C=gddot/(f^2 g)+fddot/(2f^3)-3 fdot^2/(4f^4).

No boundary term is silently used to alter the original physical
current or the selected Cauchy covariance.

For T, take f_T=N/a. Then C_T=1 and the equation is exactly
v_T''+(k^2+a^2m^2 bm)v_T=0.

For L, take f_L=N sqrt(bm/am)/a and put

    A=a(am bm)^(1/4), U=a^2m^2 am.
    C_L=m A k/sqrt(k^2+U).

The constant m does not affect C_L''/C_L. Writing b=A'/A,

    (log(C_L))'=b-U'/[2(k^2+U)]

gives the exact potential

    V_k=V_infinity+R_k,
    V_infinity=U-A''/A,
    R_k=(bU'+U''/2)/(k^2+U)-3U'^2/[4(k^2+U)^2].

A separately differentiated pump verifies the same identity.
This is not an equality to a single momentum-independent scalar
theory: the finite-k remainder is still present. Off the original
clock, sigma_T and sigma_L are generally different coordinates.
They do not change the vector fronts measured in the original
physical metric.

## Original-clock control

On the original clock N=am=bm=1, sigma_T=sigma_L, d/dsigma=a d/du,
a=(1+u^2)^2, and H=4u/(1+u^2). With z=k^2/(k^2+a^2m^2),

    V_k=a^2[m^2-z Hdot+z(1-3z)H^2],
    V_infinity=a^2[m^2-Hdot-2H^2],
    R_k=a^2(1-z)[Hdot+(2+3z)H^2].

These reproduce the frozen cosmic-time T/L oscillator potentials
after their exact time change. On |u|<=1/2, Hdot>0, so R_k>=0.
The potential extends continuously to k=0, where V_k=a^2m^2;
the canonical map C_L itself is only asserted for k>0.

Write r=1+u^2. Direct calculation gives

    U=m^2r^4, U'=8m^2u r^5,
    U''=8m^2r^6(1+11u^2),
    b=4u r, A''/A=4r^2(1+7u^2),
    bU'+U''/2=4m^2r^6(1+19u^2).

Since U>=m^2, the negative square in R_k gives the upper bound

    R_k <=(359375/4096) m^2/(k^2+m^2)
          <88 m^2/(k^2+m^2).

The numerator maximum is its endpoint value, because it is an even
polynomial with nonnegative coefficients. Similarly A''/A<=275/16,
so V_infinity>=m^2-275/16. This is positive for m=1000 and every
point of I. It is not a general positive-potential theorem for
arbitrary finite off-clock backgrounds.
