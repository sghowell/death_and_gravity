# Exact rational actual-parameter bounds

On |epsilon|=1/16, the same Gamma estimates as sector.md
give, for D(x)=xM+(1-x)^2 in [1,M],

    |T_D(1)|<70 mF^(1/8)/Q,
    |T_D(M)|<70 mF^(1/8)M^(17/16)/Q,
    |B_HL,D(-1)|<64 mF^(1/8)M^(1/16)/Q.

Use exact rational power caps mF<=U^8, 3M<=A^8,
M<=B^16. The four scalar finite-part allowances are

    quartic/OS local: 4900 L U^2/(8Q_lower^2),
    cubic sunset:    1200 g U^2 A^9/Q_lower^2,
    mixed OS mass:   2240 g U^2 B/Q_lower^2,
    heavy MS mass:   280 g U B^17/Q_lower^2.

The explicit 1/epsilon in the last term supplies the
factor 16. These are bounds on the whole meromorphic
products, not products of separately extracted finite
parts. The Cauchy-average inequality justifies adding
their absolute allowances even if cancellations would
sharpen the result.

At the actual boundary take U=A=10^25, B=10^13 and
Q_lower=144. All power-cap conditions are checked
exactly. The complete scalar allowance is approximately
8.62335717236554e265, below 1e275. It is intentionally
loose and suffices for the complete reference.

The scale logarithm log(mF^2) has the inherited exact
upper bound ell_upper<1000. Both scalar first-loop
vacuum determinants have magnitude at most

    E_scalar1=(1+M^2)(ell_upper+3/2)/(4Q_lower).

With Q_upper=256 the complete V1 is at least
63mF^4/256-E_scalar1>0. Adding the two already bounded
fermionic vacuum families to the new scalar allowance
gives, approximately,

    |V2| upper           7.77166076157809e593,
    |V2|/V1_lower upper  3.15800818248253e-206.

Strictly, these are below 1e595 and 1e-203 respectively.
The finite first source is bounded by
G(ell_upper+1)/(2Q_lower)<0.001. Its full regulated
source-square finite piece is below 7.82e-205 before
the exact triple cancellation; pi^2<16 suffices.

The large dimensionful vacuum reference and its small
relative perturbative ratio are different statements.
Neither bounds a physical omitted higher-loop term
or proves finite-gravity backreaction control.
