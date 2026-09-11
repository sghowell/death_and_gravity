# Independently fixed quadratic Newton reference

Retain the S6.168 vacuum and saturated pole-mass reference and the
S6.169 flat curvature improvement. The separately named curved
GY14-SAT8-MR/EC-N0 extension adds -m^2 R/(6Q) per Dirac copy at mu=m.
It sets this free quadratic sector's zero-field soft-metric Newton
contribution to zero; the coefficient is not adjusted to the energy
estimate. It is not the full interacting physical Newton dictionary.

Define ell=log(M^2/m^2), Dfun=M^2(ell-1)+m^2. The referenced action is
L2=-(partial M)^2(ell+2/3)/Q-Dfun R/(6Q), so

    rho2=[-Mdot^2(ell+2/3)-Dfun H^2-H Dfun_dot]/Q,
    P2=[-Mdot^2(ell+2/3)+Dfun(2Hdot+3H^2)/3
         +Dfun_ddot/3+2H Dfun_dot/3]/Q.

Since Dfun(m)=Dfun'(m)=0, Dfun''=2ell+4 and |ell|<=ellmax
=2Delta/(m-Delta)<.01 on the actual domain, Taylor's integral
remainder gives |Dfun|<=3Delta^2. Exact differentiation gives

    |Dfun_dot|<=2(m+Delta)Delta ellmax/tau,
    |Dfun_ddot|<=[2ellmax(Delta^2+9(m+Delta)Delta)+4Delta^2]/tau^2.

The full real profile has |Mdot|<=Delta/tau,
|Mddot|<=9Delta/tau^2. Inserting |H|<=2, |Hdot|<=4 gives the
uniform bounds in reference.py, with Q>144 and six active copies:

    |rho2|<1e594,  |P2|<1e595.

No inequality is differentiated. For each of36 constant copies
M=m the referenced local two-derivative tensor is exactly zero.
Their higher geometric/state remainder and Euler term are still kept.
At the actual clock bounce, the active two-derivative rho and P are
nonzero and opposite; that cancellation alone does not fix the full
tensor, scalar equation or background response.
