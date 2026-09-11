# Full sixth-order Cauchy-Leibniz coefficients

S186 bounds the actual full phase-stripped reference pair
amplitude by C sqrt(nu mu), C=1e9, on complex-time discs
of radius r=1/20000. It includes both constant actual-alpha
factors and the constrained ten-field readout. The inverse
sum phase is analytic with |g|<=G=2/(nu+mu).

Banach-valued Cauchy estimates bound amplitude and inverse
phase derivatives by j! C sqrt(nu mu)/r^j and j! G/r^j.
For arbitrary smooth test h, all product-rule terms give

    ||L^n(a h)|| <= C sqrt(nu mu) G^n
                   sum_(j=0)^n c[n,j] r^(-(n-j))||h^(j)||.

The coefficients satisfy c[0,0]=1 and

    c[n+1,j]=(2n+2-j)c[n,j]+c[n,j-1],

with missing entries zero. To see the first factor, multiply
a term of L^n(a h) by g. There are n+1 inverse-phase factors
and one amplitude factor; their total existing derivative
order is n-j. Differentiating any of these contributes the
sum of its factorial ratios, (n-j)+(n+2)=2n+2-j.
Differentiating h instead gives the preceding coefficient.

Independent direct symbolic product rules check every
coefficient for n0..6. The third row is(48,33,9,1),
the fifth is(3840,2895,975,185,20,1), and the sixth is

    (46080,35685,12645,2640,345,27,1).

The fifth and sixth sums are7916 and97423. Since r<1,
r^(-(n-j))<=r^-n. No derivative of the rapidly oscillating
actual-state mixing is used: this calculation concerns
only the analytic reference already separated in S197.
