# Unequal-mass sunset and the finite Laurent coefficient

Let e=epsilon and Q=16pi^2. Gaussian momentum integration
and a Schwinger scale integral give

    S_D=Q^-2 exp(2 gamma_E e) mF^(4e) Gamma(2e-1) A(e),
    A(e)=integral_simplex U^(e-2) X^(1-2e),
    U=xy+xz+yz, X=a x+b y+c z.

Here a,b,c are positive squared masses. The momentum
representation and all massive tadpoles/bubbles have a
common convergent dimensional strip 1<Re(e)<2 before
meromorphic continuation. The standard scaleless constant
traces are removed by their dimensional identity.

Divide the simplex into its six ordered coordinate sectors.
For each labeled mass permutation set

    (x,y,z)=(1,r,rt)/(1+r+rt),  0<r,t<1,
    V=1+t+rt, B=a+r b+r t c.

The Jacobian is r/(1+r+rt)^3. All powers of 1+r+rt cancel,
leaving r^(e-1) V^(e-2) B^(1-2e). Subtract its r=0 value:

    A_sector(e)=a^(1-2e) J(e)/e
       +integral dr dt r^(e-1)
          [V^(e-2)B^(1-2e)-(1+t)^(e-2)a^(1-2e)],
    J(e)=(2^(e-1)-1)/(e-1).

The difference is O(r), uniformly on compact e sets.
The subtracted integral is therefore analytic for
Re(e)>-1. This explicitly isolates every corner pole;
the Gamma factor retains the overall pole. Both are
kept until the full finite Laurent coefficient.

For the actual scalar masses 1,1,M with M>=2, use
rho=1/16 and |e|=rho. All positive mass forms have
unambiguous real logarithms for their complex powers.
For F=V^(e-2)B^(1-2e),

    |partial_r F| <=
      [(2+rho)+(2/3)(1+2rho)] (3M)^(1+2rho)
      <3(3M)^(9/8).

The six compact remainders are bounded by
(96/5)(3M)^(9/8), and the six corner references by
96 M^(9/8). Hence |A(e)|<120(3M)^(9/8).

For 7/8<=Re(z)<=9/8, the Gamma integral gives
|Gamma(z)| <= 8/7+2/exp(1)<2. The Gamma recurrence and
0<gamma_E<1, together with exp(1/8)<2, then imply

    |exp(2 gamma_E e) Gamma(2e-1)|<40,
    |exp(gamma_E e) Gamma(e-1)|<70,
    |exp(gamma_E e) Gamma(e)|<64.

The explicit denominator lower bound for the first is
|2e(2e-1)|>=(1/8)(7/8)=7/64. The elementary estimates
use exp(1)>5/2 and exp(x)<=1/(1-x) for 0<x<1;
no numerical special-function inequality is assumed.

Combining the sector and Gamma estimates gives

    |S_D|<4800 mF^(1/4)(3M)^(9/8)/Q^2.

The finite Laurent coefficient is the circular average
of the meromorphic expression. Its magnitude is at most
this supremum; pure poles average to zero. This bounds
the complete finite coefficient, including every
positive-epsilon numerator times a pole.
