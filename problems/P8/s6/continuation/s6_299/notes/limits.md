# Ordered exponentially fine detector limit and explicit wrong-order test

First remove e at fixed positive x=E/nu and fixed kappa. Write
a=A/kappa,Delta=D/kappa with A>=0 and D fixed by the hard kinematics.
The complete leading factor is exp(D/kappa)*F(A/kappa)*x^(A/kappa).

For x=exp(-kappa*chi),chi>=0, this becomes
exp(D/kappa)*F(A/kappa)*exp(-A*chi).
Taking kappa->infinity yields exp(-A*chi). The power is kept intact.
A first-order truncation1-A*chi can even be negative; it is not a rate
approximation uniformly in this regime. For A>0 the iterated limits
x->0 then kappa->infinity give0, while reversing them gives1.
These are statements about the same leading-soft factor, not a choice
to alter the physical prescription.

An explicit regulator-order check is available. Instead take
e=rho/kappa withrho>0 fixed at the same time asx=exp(-kappa*chi).
For the auxiliary Poisson representation,

lambda -> A*exp(-2rho*chi)/(2rho),
2e*N ->0 in L2,
lambda-a/(2e) -> A*(exp(-2rho*chi)-1)/(2rho).

The resulting soft factor is
exp[A*(exp(-2rho*chi)-1)/(2rho)],
not exp(-A*chi). Its logarithmic excess over the ordered result is

A*integral_0^chi [1-exp(-2rho*u)]du,

which lies in[0,A*rho*chi^2]. Thusrho->0 recovers the ordered limit;
arbitrary finite rho does not. This path calculation is not a replacement
for the mandated fixed-resolution regulator removal. It exposes why a
regulator-small assertion without a resolution-uniform estimate is unsafe.

The original finite conversion error bound remains uniform after e has
been removed. S296's finite-e bound itself contains lnE; neither that
bound nor this packet silently asserts a uniform arbitrary joint limit.
No physical hard kinematic parameter, reference scale or detector value
is selected to force a P8 verdict.
