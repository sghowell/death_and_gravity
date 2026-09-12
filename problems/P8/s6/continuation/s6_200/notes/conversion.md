# Exact cubic division and spatially nonlocal equal-time terms

For fixed p set b=s+q and expand the complete numerator
N_i(w-q,p)=sum_(k=0)^3 n_ik w^k. For any cubic polynomial N,
ordinary division gives

    N(w)/(b-w)
      =sum_(r=0)^2 w^r sum_(k=0)^r n_k/b^(r-k+1)
       +w^3 N(b)/(b^3(b-w)).

Multiplication by rho_i(s)/s^3 and integration is legitimate:
each displayed conversion coefficient has an absolute s^-2
or better ultraviolet majorant. The residue is
N_i(b-q,p)=N_i(s,p)=s Q_i(s,p), yielding precisely the B6
kernel in FORMULATION.md. No component or endpoint is discarded.

Equivalently A_r=partial_w^r F(w-q,p)|_(w=0)/r!.
The expansion center is z=-q, not z=0. The independent numerical
interpolation initially confused these centers; direct frequency
derivative tests now check all four centered coefficients for both
spin sectors at four spatial transfers. Tests retain the original
tight integral comparisons.

For a unit-norm spatial TT direction transverse to p, the tensor
numerator reduces to z^3 and the scalar channel vanishes. In this
channel

    A0=-q^3 integral rho2(s)/(s^3(s+q)) ds.

For |q|<4m^2 its coefficient of q^(3+n) is
(-1)^(n+1) integral rho2(s)/s^(4+n) ds. Every such moment is
strictly positive and finite, so infinitely many Taylor coefficients
are nonzero. Therefore the equal-time conversion is not a finite
spatial differential counterterm. Setting it to zero or absorbing it
into arbitrary local constants would change the full spatial kernel.

The two sides are explicitly defined continuum nonlocal representatives.
This calculation neither selects the physical local matching
polynomial nor identifies the original one-leg contact regulator with
the two-leg memory regulator.
