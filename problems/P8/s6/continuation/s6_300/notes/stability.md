# Uniform multiplicity bound and correct nested IR pairing

For mu1 put c=2x(1-x) in[0,1/2] and d in[1,7].
Differentiating the Feynman representation gives an integrand numerator

2d+c*(d^2-2d+1/2)
=2(d-1)+c(d-1)^2+7/4+(1/4-c/2)>0.

The negative term in the original derivative also gives
F'(d)<=2d*integral1/[1+c(d-1)]<=14.
Together with the six-dot shift bound, the massive contribution to
K changes by at most2*14*48R=1344R.

On a_i in[1/4,4], |a_i ln(2a_i)|<9:
ln2<3/4 follows from exp(3/4)>1+3/4+(3/4)^2/2>2.
The complete mixed contribution is bounded by2*4*9R=72R.
On b in[0,2], |b ln(2b)|<=3, including its negative minimum
-1/(2exp1) and its positive endpoint2ln4<3.
The radiation/radiation term is at most3R^2.

Consequently, for0<R<=1/8,

|K_sigma-K_elastic|<=1416R+3R^2<1440R,
|a_sigma-a_elastic|<40R/kappa,

using a=K/(4pi^2*kappa) and pi>3. At R0 the difference is exactly zero.
The estimates have no N dependence and survive collinear splitting
and the finite-energy-measure limit.

For a fixed radiative state the next infinitesimal graviton contributes
a real rate pole a_sigma/(2e). It must pair with virtual rate pole
-a_sigma/(2e). Using the elastic virtual term instead leaves
(a_sigma-a_elastic)/(2e). A tiny bound on its numerator is NOT a finite
remainder at fixed R when e->0. The universal pole is known here; the
finite radiative hard amplitude and Delta_sigma are not.

After state-correct pairing consider only the KNOWN factors
x^a*F(a), F(a)=exp(-EulerGamma*a)/Gamma(1+a). Assume nu1 and
R<=x<=1/8. Then
|delta ln(x^a)|<=40R|lnx|/kappa<12/kappa,
because -xlnx is increasing up to1/8 and ln8<9/4.
The old elastic index is<4/(5kappa), and the new index is<6/kappa.
The Weierstrass derivative
|dlnF/da|=sum_j a/[j(j+a)]<=pi^2*a/6<2a
therefore gives |delta lnF|<60/kappa^2.

At original kappa=10^800 the combined logarithmic change is<13/kappa.
Since this is<1/2, the multiplicative ratio differs from one by
<26/kappa<3*10^-799. The resolution power is unexpanded. This is a
bound on those two known factors, NOT a full detector-rate error.
The unknown finite angular and hard radiative terms remain explicit.
