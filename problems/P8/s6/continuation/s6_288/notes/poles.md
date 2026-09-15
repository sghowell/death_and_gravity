# Complete dimension-dependent nonbox physical principal parts

Write a+b+c=4mu, D=4+2EP, and
VD(a)=a^2-4mu*a+4mu^2(D-3)/(D-2).
The complete tree channel is
N_D(a,b,c)=4mu^2(D-3)/(D-2)-2mu*a-bc,
and A_tree,D=-sum N_D/(kappa*a). In particular N_D(0,b,4mu-b)=VD(b).

Let B0 denote the raw massive zero-momentum bubble:
B0=-Gamma(1-EP)(4pi*nu^2)^(-EP)*mu^EP/EP.
Its complete parameter integrals give
C0mumu(0)=B0/(2mu),
Bmm_prime(0)=-EP*B0/(6mu), A0(mu)=mu*B0/(1+EP).
The first equality is an equality in the common dimensional variable:
the C0 pole is infrared and the B0 pole ultraviolet. Their physical
origins are not identified.

For the on-shell mixed bubble,
Bmix(mu)=B0/(1+2EP) and
Bmix_prime(mu)=-B0/[2mu(1+2EP)].
Substitution into the entire S283 stress self-energy gives
Sigma(mu)=0 and
Sigma_prime=mu*B0*r(D)/(16pi^2*kappa),
r(D)=2(D-1)/[(D-2)(D-3)].
Its finite expansion is exactly
mu[-3/EP+7+3log(4pi*nu^2/mu)-3Gamma_E]/(16pi^2*kappa).
All four external legs give-2Sigma_prime*A_tree,D.
The tree must remain D dependent when multiplied by its poles.

Let b2,b1 be the a^-2,a^-1 coefficients of the exact S284 Bmm
coefficient, and c1 the a^-1 coefficient of its C0mumu coefficient.
The entire nonbox-plus-LSZ principal part is

 B0 { b2/a^2+
 [b1-EP*b2/(6mu)+c1/(2mu)+2mu*r(D)*VD(b)]/a }.

Both C00mu and B00 coefficients are regular at a0. The entire crossed
S284 Gram choice is regular there as well. The code compares the
literal full-D frozen coefficients, not only their D4 limits.

Define the crossing basis
E3=mu^3 sum1/a, E4=mu^4 sum1/a^2,
E5=mu sumbc/a, E6=mu^2 sum(b-c)^2/a^2.
The complete principal parts in all three channels are B0*R_D,
where R_D=sum rj Ej and

 r3=-2(2D^3-78D^2-5D+156)/[3(D-2)^2(D-1)(D+1)],
 r4=-8(2D-7)/[(D-2)^2(D-1)(D+1)],
 r5=(D^4-41D^3+86D^2+56D-108)/
     [3(D-3)(D-2)^2(D-1)(D+1)],
 r6=1/[(D-1)(D+1)].

At D4 these are [164/15,-2/15,-73/15,1/15].
Their EP-linear coefficients are
[-4879/225,-28/225,3128/225,-16/225].
Thus if B0=-1/EP+L+O(EP),
the finite principal part is L*R4-R_EPprime. Omitting either the
evanescent coefficients or the evanescent tree changes this answer.

The infrared C0 principal part is-2mu*V(b)/(a*EP_IR); the four
external-leg IR part is+2mu*V(b)/(a*EP_IR). They cancel. The remaining
meromorphic pole is ultraviolet. This does not cancel the entire
finite-transfer soft factor or the box Coulomb phase.

The S284 Gram-removed bubble UV polynomial plus the full LSZ UV
contribution8mu*sum N4/a has crossing coefficients

 [203/40,-169/3,164/15,-2/15,-73/15,1/15]

in [sum a^2,mu^2,E3,E4,E5,E6]. Subtracting B0*R_D removes all four
meromorphic coefficients. Its chosen cut representative's remaining
UV coefficient is203/40 sum a^2-169/3 mu^2. Unknown regular tadpole
and counterfunctional terms may alter the local matching description;
this is not asserted to be a full-source physical beta function.
