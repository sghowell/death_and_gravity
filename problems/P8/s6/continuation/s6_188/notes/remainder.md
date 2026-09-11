# Actual finite-amplitude modewise covariance remainder

The unperturbed reference energy is
E=Z^T M0(t)Z/2. Hamiltonian evolution cancels its
generator terms, so E'=Z^T M0'(t)Z/2.
The transverse K rate is -H, the longitudinal K terms
have rates -H and -3H, and the positive V terms have
rates H and -H. Hence |E'|<=3|H|E.

On the actual slab, a falls from a0=25/16 to1 and
returns. Thus integral|H|<=2 log a0 between any two
slab times. The reference propagator satisfies

    ||S0(t,s)||M0(s)->M0(t) <= a0^3
         =15625/4096 <4.

This bound is uniform in all internal momenta,
including zero. Also
||M0^1/2 J M0^1/2||op=omega(t), because K0 V0=omega^2 I.

Take gamma(t)=epsilon Gamma(t), with Gamma symmetric
tracefree, smooth, compact after the initial neighborhood
and sup_t||Gamma(t)||op<=1. Directions at different times
need not commute. Set delta=|epsilon| and
omega_max=sqrt(m^2+|k|^2). In the fixed initial energy
metric, U=S0(t0,t)S_epsilon(t,t0) solves U'=B_epsilon U.
The two reference propagator factors and the full
exponential vertex estimate give

    ||B_epsilon(t)|| <=16 omega_max
                              [exp(delta ||Gamma(t)||op)-1].

On this unit slab define A=16 omega_max(exp(delta)-1).
Then integral||B_epsilon||<=A, ||U||<=exp A and
||U-I||<=exp A-1. Let U1=integral B1 be the exact
derivative at zero. The single-vertex Taylor tail
and the entire time-ordered Dyson tail give

    ||U-I-epsilon U1|| <= R,
    R=8 omega_max delta^2 exp(delta)+A^2 exp(A)/2.

The first term uses exp(x)-1-x<=x^2 exp(x)/2;
the second bounds every ordered term of degree at
least two by the corresponding scalar exponential.
No omitted formal-order symbol replaces the remainder.

Let C_init,E be the initial covariance in that energy
metric. It is positive and tr C_init,E<=108 nu, where
nu=omega(t0), by the actual-state bound in state.md.
Writing X=U-I and R_U=U-I-epsilon U1 gives exactly

    C_epsilon,I-C_init-epsilon C1,I
       =R_U C_init+C_init R_U^T+X C_init X^T,

Therefore

    ||covariance remainder||trace
       <=tr(C_init,E)[2R+(exp(A)-1)^2]
       <=108 nu[2R+(exp(A)-1)^2].

Returning to the instantaneous background energy metric
costs at most16, the square of the reference norm gain.
The normalization is a positive initial energy trace,
not a vanishing background-density denominator.

For the illustrative mode class omega_max<=2000 and
delta<=1e-8, take Abar=32*2000*1e-8<1/2. The elementary
exp(x)<=1/(1-x) for0<=x<1 gives the rational bound

    Rbar=8*2000*delta^2/(1-delta)
                       +Abar^2/[2(1-Abar)],
    2Rbar+[Abar/(1-Abar)]^2 <1e-6.

This is an interaction-picture covariance error relative
to the initial energy trace, not a relative error in
every component or in a renormalized stress.

The general all-momentum formula grows with omega_max.
It is NOT a uniform integrable majorant for removing a
stress-response regulator. The mode benchmark is not a
physical cutoff. UV subtraction, its finite contacts
and a continuous subtracted remainder remain required.
