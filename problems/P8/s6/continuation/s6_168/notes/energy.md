# Full momentum and full potential bounds

For both exact Hadamard states at every real time,
the subtracted-mode energy remainder has the allowance

    R_state=8NC Delta/(3pi^2 tau^4 m0),
    R_frames=24N1024^2 Delta^2/(5pi^2 tau^4 m0^2),
    R_Taylor=NDelta^4/(672pi^2 tau^4 m0^4),

where C=2^42 and m0=.99mF. The factors include two
occupied helicities and radial measure1/(2pi^2).
The required complete radial integrals are

    integral p^3/E^5 =2/(3m0),
    integral p^4/E^7 =1/(5m0^2),
    integral p^6/E^11=2/(63m0^4).

Using pi^2>9, their actual allowances are approximately

    2.36931798577993266e410,
    3.08121505968778696e401,
    8.36534214389906295e386.

Their sum is below 1e411. This is not an adiabatic
state substituted for the exact evolution.

Let u=ychi/mF. Subtract the fixed zero-field vacuum
constant and the saturated fixed pole-mass reference.
The full paired potential is

    -.5[fF(1)-fF(0)]chi^2
    -(8/3)NY^2 chi^4/Q
    +mF^4/Q sum_even_k>=6 c_k u^k,

    c_k=48N/[k(k-1)(k-2)(k-3)(k-4)]>0.

The coefficients decrease from c6=N/15. For
|ychi|<=Delta the positive entire tail is at most

    NDelta^6/[15Q mF^2(1-(Delta/mF)^2)].

This is the full convergent potential, not just its
quartic Taylor polynomial. No sign cancellation is
used between the remaining terms. The absolute
mass-anchor and quartic allowances are

    2NDelta^2/(3Q), 8NDelta^4/(3Q).

The first expression includes the fixed m_Phi^2=1
reference factor; all energy terms have dimension four.

For M=mF +/- Delta s(t/tau), |Mdot|<=Delta/tau and

    |log(M^2/mF^2)|<=2(Delta/mF)/(1-Delta/mF).

The actual logarithm bound is below1/100. With the
2/3 dimensional finite term, the absolute restored
derivative energy is at most NDelta^2/(Q tau^2).
The same factor-below-one bound holds throughout the
accepted Delta/(.99mF)<1/100 parameter domain.

Use mF=10^200, Delta<=3e197, tau=10^-100, N=6,
and Q>144. The four local allowances are approximately

    mass anchor: 2.5e393,
    quartic: 9e788,
    entire higher potential: 2.02501822516402648e782,
    derivative finite part: 3.75e593.

Adding every local allowance and every state remainder
gives 9.00000202501822516402648e788, strictly below
1e789 and below1e-11 times kappa*m_Phi^4.

This is an absolute allowance for the explicit
one-loop quadratic/reference contribution, not a
relative error against the zero bounce density,
not full stress, and not the interacting parent.
