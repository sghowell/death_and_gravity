# Literal external self-energy and complete soft pole

For k=p-r the exact off-shell Ward identity is
k:T=(p^2-mu)r-(r^2-mu)p. Direct tensor contraction in D gives

    T:P:T=2p^2r^2+4mu p.r-2Dmu^2/(D-2)
      =2(p^2+mu)(r^2-mu)+4mu p^2-4mu^2/(D-2)-2mu k^2.

The source's two vertex factors give the normalized self-energy in the
formulation. The massless tadpole and scalar-graviton seagull are scaleless
and vanish in the stated dimensional convention. No massive tadpole is
dropped.

Let C=Gamma(eps)/rGamma*(nu^2/mu)^eps. At the exact tree pole,

    B0(mu;0,mu)=C/(1-2eps), A0(mu)=mu*C/(1-eps),
    4mu^2-4mu^2/(D-2)=2mu^2(1-2eps)/(1-eps).

Their substitution proves Sigma(mu)=0 EXACTLY in D. Dropping the
evanescent trace before the loop would leave an incorrect finite result.
This is a statement about the specified pure-gravity graph sector.

Differentiating the full bubble parameter integral at the pole gives

    eps*integral_0^1 x^(-1-2eps)(1-x)dx=-1/[2(1-2eps)],
    B0'(mu)=-B0(mu)/(2mu).

The exact derivative is
mu(3-2eps)B0(mu)/(1-eps), with prefactor1/(16pi^2 kappa).
Its expansion, including the raw Gamma_E/log4pi conversion, is recorded
in the executable packet. The4mu/eps_UV term comes from4mu B0;
the-mu/eps_IR term comes from the derivative master. They cannot be
interpreted as the same physical divergence merely because analytic
dimensional notation combines them.

The inverse convention is p^2-mu+Sigma, so Z=1/(1+Sigma').
Four external scalar legs multiply the proper amplitude by Z^2:
the first-loop correction is-2Sigma'*A_tree. Its IR contribution for
EP=-eps_IR>0 is-mu*A_tree/(8pi^2 kappa EP).

In the common four-point master basis, the pairwise box/C0mumu pole
has coefficient4kappa A_tree*sum V(a)M(a), inside
1/(16pi^2 kappa^2 EP). Adding the four-leg coefficient
-2mu*kappa*A_tree yields

    A_tree/[8pi^2 kappa EP] * [2sum V(a)M(a)-mu],

the ENTIRE accepted S278 soft pole. This includes the self term that
four-point normal cuts alone do not determine.

The raw finite residue is not a complete renormalized massless-pole anchor.
Proper vertex terms, other sectors, vacuum/gravitational counterterms and
finite physical matching are still required; none is set to zero here.
