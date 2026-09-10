# Entire scalar finite remainder and its nonzero pole

The compact remainder in scalar.md is

    D3(e)=Integral_0^1 z^(2e-3)(1-z)^(3/2-e) h3(z/T) dz,
    h3(k)=k^3/(1-k).

It is holomorphic for |e|<1/8. Near zero, h3(z/T) removes
all three negative powers of z. Near one, the exponent is
at least 11/8. Regulator derivatives add integrable logarithms.
This uniform domination permits coefficientwise extraction.

The prefactor is Gamma(e-1)H(e)=-1/e+4log2-3+O(e).
D3(0) is nonzero. Thus the overall simple pole
-T^2 D3(0)/e must be subtracted, while its finite prefactor
product remains. The exact finite remainder is

    delta_V=T^2 Integral_0^1 (1-z)^(3/2)/z^3 h3(z/T)
       [4log2-3-2logz+log(1-z)] dz.

For T>=16, h3(k)<=2k^3 follows from
2k^3-h3=k^3(1-2k)/(1-k)>=0. Also
|4log2-3|<=3 and both unit-interval negative-log moments equal
one. Dropping the harmless factor (1-z)^(3/2)<=1 gives

    |delta_V| <= T^2*(2/T^3)*(3+2+1)=12/T.

No finite mass-ratio term is rounded to zero. The exact
rational enclosure uses pi<4, so
13/4+pi^2/8<21/4 and 47/18+pi^2/12<4.
Together with the absolute massless reference,

    |V_scalar,paired| <= (6Y_upper/Q_lower^2)
                         [19m^4+21T/4+4+12/T].

Q_lower<=144 is a conservative lower bound on physical 16pi^2.
The gauge enclosure is
42 a_upper (4/3) 18m^4/Q_lower^2.
All input bounds are explicit and exact rational values.
