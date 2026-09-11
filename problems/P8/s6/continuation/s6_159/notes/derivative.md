# Two-light-mass sunset derivative without a full heavy-mass loss

Differentiate the S6.158 corner-subtracted sunset with
respect to both light squared masses a at a=1. In the
three heavy positions, each occurring twice among six
labeled sectors,

    B=M+a r+a rt,  C=r+rt;
    B=a+M r+a rt,  C=1+rt;
    B=a+a r+M rt,  C=1+r.

Here C=partial_a B, 0<=C<=2 and 0<=partial_r C<=2.
With V=1+t+rt the differentiated integrand is

    Fa=(1-2epsilon) C V^(epsilon-2) B^(-2epsilon).

The corner value is differentiated as well. All masses
are positive; differentiation first holds in the common
convergent strip and extends through the same compact
corner subtraction. The finite Laurent coefficient is
taken only after this operation, at fixed scale.

On |epsilon|=rho=1/16, at a=1 and M>=2,

    |Fa| <= (9/4)(3M)^(1/8),
    |partial_r Fa| <
       (9/8)[2+2(2+rho)+8rho M](3M)^(1/8)
       <10 M(3M)^(1/8).

Consequently the corner difference is at most
min(10Mr,5)(3M)^(1/8). Bounding its derivative uniformly
by a full power of M over the whole interval would be
unnecessarily weak. Instead split r at r0=1/(3M):

    small-r multiplier: (10/3)/(1-rho)=32/9,
    large-r multiplier: 5/rho=80.

Both multiply (3M)^(rho+1/8)=(3M)^(3/16).
Their sum is below 84 per sector, hence below 504 for
six sectors. The six differentiated corner references
are bounded by 108(3M)^(1/8). The total is below
700(3M)^(3/16).

The inherited Gamma/exponential factor is below 40.
Thus

    |partial_a S_D(a,a,M)|_(a=1)
      <28000 mF^(1/4)(3M)^(3/16)/Q^2.

The finite Laurent coefficient obeys the same
Cauchy-circle bound. This controls the entire regulated
mass derivative and its pole products, not a finite
part differentiated after changing its reference.
