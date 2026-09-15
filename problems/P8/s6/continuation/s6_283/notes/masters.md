# Complete scalar master normalization and uniform box reduction

The named integrals use D=4-2eps, measure
nu^(2eps)/(i*pi^(D/2)*rGamma), and propagators with+i0.
Here rGamma=Gamma(1-eps)^2 Gamma(1+eps)/Gamma(1-2eps).
The explicit code supplies the Feynman-i0 limits, not real-part-only formulas.

For the alternating box with internal masses0,mu,0,mu and four external
mass squares mu, the exact Symanzik form is

    F=mu(alpha2+alpha4)^2-s alpha1 alpha3-t alpha2 alpha4.

Let b=alpha2+alpha4, alpha2=bx,alpha4=b(1-x),
alpha1=(1-b)y,alpha3=(1-b)(1-y), and r=b/(1-b).
The Jacobian is b(1-b). With A=mu-tx(1-x), B=-sy(1-y),
the entire parameter integrand becomes

    r(1+r)^(2eps)/(A r^2+B)^(2+eps).

First work at s,t<0 and Re(eps)<0, where all integrals converge.
Replacing(1+r)^(2eps) by1 integrates r exactly to
B^(-1-eps)/[2A(1+eps)]. The y beta integral and checked gamma recurrence
then give4M(t)/(s eps)*(nu^2/(-s))^eps.

This replacement loses no finite term. For |eps|<1/8,
|(1+r)^(2eps)-1| is bounded by
2|eps|log(1+r)(1+r)^(1/4). A is bounded above and below by positive
constants. For r<=1 and y near0 or1, B is comparable to the distance
to that endpoint. Integrating y in the difference leaves a bound
C|eps|r^(-1/4), integrable in r. Away from those endpoints compact
domination applies. For r>=1, the difference is bounded by
C|eps|r^(-5/2)log(1+r), uniformly integrable. Thus the entire difference
is O(eps), locally uniformly in the spacelike invariants. Analytic
continuation with the Feynman prescription establishes

    I4(s,t)=4M(t)/s[1/eps+log(nu^2/(-s))]+O(eps).

This agrees with Box14 of Ellis and Zanderighi,
https://arxiv.org/html/0712.1851 , but the normalization and absence of
an extra finite term follow from the explicit derivation above.

For C0mumu, group its two massive parameters into b. Its denominator
is b^2[mu-sx(1-x)-i0], and the measure is b db dx. The b integral
is-1/(2eps), including the triangle's overall minus sign:

    C0mumu=Gamma(1+eps)/(2eps*rGamma)*nu^(2eps)
                 integral_0^1[mu-sx(1-x)-i0]^(-1-eps)dx.

Since Gamma(1+eps)/rGamma=1+O(eps^2), this supplies its complete pole
and finite logarithmic parameter integral. C00mu is finite at nonzero
spacelike s: its denominator is mu*b^2-s(1-b)^2x(1-x), with measure
(1-b)dbdx. The corner bound is the integrable
integral dbdx/(b^2+x), and similarly at x=1. Both bubbles have the
ordinary one-parameter logarithmic representation.

Convert eps=-EP to the retained raw convention before comparing cuts:
rGamma(-EP)(4pi)^(-EP)=1+EP(Gamma_E-log4pi)+O(EP^2).
This fixes every displayed Gamma_E and4pi term. It is not an implicit
MSbar finite counterterm choice.
