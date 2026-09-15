# Complete compact-window Laurent identity and remainder

Let mu>0,n>4mu,kappa>0,g!=0,0<L<=1/2,0<e<=e*=1/8 and W a fixed real C1
test weight on[n,n(1+L)]. Set y=(s-n)/n and
H_e(y)=n^(2e)K_e(n(1+y))W(n(1+y)). Then
int_n^{n(1+L)}rho_e(s)W(s)ds=int0^L y^(-1+2e)H_e(y)dy.
The exact subtract-and-add identity is
I_e=H_e(0)L^(2e)/(2e)+int0^L y^(-1+2e)[H_e(y)-H_e(0)]dy.

Writing H0=H_e|0,H1=d_e H_e|0 gives
I_e=H0(0)/(2e)+H0(0)lnL+H1(0)/2+
int0^L[H0(y)-H0(0)]/y dy+R_e.
Here H0(0)=K0(n)W(n) and
H1(0)=[K1(n)+2ln(n)K0(n)]W(n).
The logs therefore combine as n/(16pi nu^2), a dimensionless ratio.
Using a dimensionful windownL withK1 gives exactly the same finite term.

Define the following suprema on0<=e<=1/8,0<=y<=L:
M2=sup|d_e^2[H_e(0)L^(2e)]|,
M10=sup|d_y H_e(y)|,
M11=sup|d_e d_y H_e(y)|.
Then

abs(R_e)<=e*[M2/4+2M10 L(1-lnL)+M11 L].

For the first term Taylor's integral remainder divided by2e givese M2/4.
For the regular term apply the mean-value theorem in e, then in y to
both H_e(y)-H_e(0) and its e derivative. Use y^(2e)<=1 and
int0^L|lny|dy=L(1-lnL). This proves the complete displayed bound.

All suprema are finite: throughout the compact massive domain,
1-B(s)x^2>=4mu/[n(1+L)]>0. Angular derivatives through second e order
and one s order are dominated by finite multiples of
(1-x^2)^2(1+abs(ln(1-x^2))^2), which is integrable. Gamma factors are
smooth on[0,1/8], and W,W' are bounded. These are explicit conditional
supremum bounds; they are not evaluated physical detector constants.

The complete kernel lies between(1-x^2)^2 and1, so8/15<=F0<=1 and
2ln2-2<=L_B<=0, in particular-2<L_B<0 for the actual massive interior.
At mu1,nu^2=1,4<n<10^200 usepi<4,0<EulerGamma<1 andln10<5/2.
The gamma bound follows from harmonic sum/integral inequalities;
the ln10 bound already follows from the positive exponential sum
throughdegree6 at5/2 exceeding10. Thusln(16pi n)<202ln10<505 and
abs(K1(n))/K0(n)<509*15/8<1000.
This bounds the finite evanescent coefficient only, not the IR pole.

K0(n)=g^2(n-4mu)^2F0(1-4mu/n)/(8pi kappa n)>0.
If W(n)>0 the unpaired Laurent pole is positive, and the naive D4
integral has a logarithmic lower-endpoint divergence. A smooth positive
dispersion weight is an example, but this single unpaired cut does not
define the physical dispersion observable.

H is unstable in the original theory. Virtual/pole terms, heavy width
and physical IR limits must be paired before claiming a finite result.
The window identity is bookkeeping, not an arbitrary finite subtraction,
new counterterm, physical-divergence verdict, or a P8 no-go.
