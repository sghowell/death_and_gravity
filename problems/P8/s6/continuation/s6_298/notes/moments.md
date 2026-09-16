# Complete positive moment, finite tail and actual hierarchy

Set F(z)=mu(1-z)^2+nz. The complete generated derivative is

f1'(0) = g^2/(96pi^2) integral_0^1 z^2(1-z)^2/F(z)^2 dz.

The separate light and heavy moments have numerators z^2(1-z)^3 and
z^3(1-z)^2. Their sum is the full two-triangle expression, not a selected
threshold estimate. The positive Stieltjes representation inherited from
S286/S290 has support4mu and4n. Its unsubtracted moment is

Pi'(mu)=g^2/(16pi^2) integral_0^1 z(1-z)/F(z) dz.

For n>=mu>0, nz<=F<=n. Therefore

g^2/(2880pi^2 n^2) <= f1'(0) <= g^2/(288pi^2 n^2),
Pi'(mu) <= g^2/(32pi^2 n).

The strict integrand inequalities apply away from endpoints. The known
formal-graph tail satisfies
integral_U^infinity Im f1(T)/(pi T^2)dT <=Pi'(mu)/U,
by positivity and1/T^2<=1/(UT). This extends the known loop graph,
not a claim that it is the unknown full UV measure.

For r=b/a<1 the full positive even-P2 series obeys

f2 = g^2/(4pi a) sum_(k>=1) k r^(2k)/[(2k+1)(2k+3)]
   <=g^2 b^2/[60pi a(a^2-b^2)].

The coefficient inequality follows from
1/15-k/[(2k+1)(2k+3)]
=(4k-3)(k-1)/[15(2k+1)(2k+3)]>=0.
For the light cut this bounds
2 Im f1/(pi T^2) <=g^2(T-4mu)^2/(240pi^2 n^3 T^2).
Thus if4mu<U<=4n, its partial slope is
<g^2(U-4mu)/(240pi^2 n^3), a fraction<6U/n of the full slope.

For the actual hierarchy put a=mu/n and delta=10^-98, so a<delta^2.
Let J be n^2 times the total parameter integral. Rewriting its integrand
as (1-z)^2/[1+a(1-z)^2/z]^2 shows

1/3-7delta/5 < J <1/3.

Indeed1-(1+y)^-2<=2y for y>=0. On0<z<delta the missing integral is
<=delta; on the remainder it is <=2a/(5delta).
For H, the scaled heavy integral, the same pointwise inequality yields

1/12-2a/5 < H <1/12.

It follows algebraically that |H/J-1/4|<2delta for delta<1/10.
The complete f1' equals g^2/(288pi^2 n^2) with relative error<10^-97.
The original4<T<=16 window captures less than96/n<10^-195 of the slope.
Both bounds are exact; 280-digit quadrature is an independent calibration.
The equal-mass inverse-light-mass enhancement is not this source's formula.
