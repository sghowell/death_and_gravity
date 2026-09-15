# Whole-topology spacelike continuity without Gram division

The proof concerns the pole and finite Laurent coefficients of the
specified one-loop ordinary F1 after covariant local UV subtraction.
It does not prove a finite derivative, an unregulated finite-transfer
observable or an exchange of the transfer and IR limits.

Take mu=nu=1, t=-tau,0<tau<=1 and0<=EP<=1/4. On-shell Feynman parameter
representations are the boundary from p^2<mu with the same i0
prescription. The spacelike parameter denominators below are positive
in the open simplex. A positive off-shell gap only strengthens the
bounds, so the stated on-shell boundary is justified.

## A direct nonsingular form-factor component

Let q=(0,sqrt(tau),0,0),
P=(sqrt(mu+1+tau/4),0,1,0), p=P-q/2, pprime=P+q/2.
Both scalar momenta have square mu and P.q=0. The conserved on-shell
tensor is
Gamma=2PP F1+(qq-eta*t)F2.
Thus Gamma02/(2P0P2)=F1 and (2P0P2)^2=4(mu+1)+tau.
This coefficient stays away from zero. A singular F2 is not divided
into F1, and no inverse transfer or tensor Gram determinant is used.

## The two-h one-Phi triangle master

Put h=(1-v^2)/4 and
A=z^2+tau(1-z)^2 h, with measure(1-z)dz dv on the unit square.
The positive scalar parameter integral is
I_EP=integral (1-z) A^(-1+EP) dz dv.
Set w=z/(1-z). Its full integrand becomes
(1+w)^(-1-2EP)(w^2+tau*h)^(-1+EP).
Dropping the first factor and rescaling w=sqrt(tau*h)*u gives

I_EP <=tau^(-1/2+EP)
       integral_0^infinity(1+u^2)^(-1+EP)du
       integral_0^1 h^(-1/2+EP)dv
     <=12 tau^(-1/2+EP).

Indeed the first integral is at most1+integral_1^infinity u^(-3/2)du=3.
The second is at most integral h^(-1/2)dv=pi<4.
Gamma(1-EP)<4/3+1/e<2 by splitting its positive defining integral at1,
and(4pi)^(-EP)<=1. The NEGATIVE raw scalar master therefore obeys
tau*abs(C00mu_raw)<=24sqrt(tau), uniformly on this EP interval.
A scalar-master bound alone is not yet the vertex result.

## Entire numerator and rank-zero terms

The retained EH GammaGamma cubic has exactly two derivatives, each
carrying one of k,k+q,q. Each literal scalar stress is at most linear in
the internal k. Therefore the entire numerator is rank at most4 and
belongs to the soft ideal(k,q)^2, component by component. This statement
uses the full EH action, not a selected helicity numerator.

After the Feynman shift k=L+z*p-y*q, the L0 coefficient belongs to(z,q)^2.
This is an algebraic ideal statement: substituting into any quadratic
soft monomial cannot produce a constant or first-degree z,q term.
The generic coefficient polynomial in the executable check verifies
the complete degree pattern. Ordinary isotropic tensor integration of
L0,L2,L4 introduces only nonsingular rational functions of D, not Gram
or transfer denominators. Odd L moments vanish.

The finitely many remaining external-momentum coefficients are smooth
and bounded in the displayed Breit family. For each L0 monomial:

- At least two q powers give tau*I_EP<=12sqrt(tau).
- One q and at least one z are bounded by
  sqrt(tau)*(17/8-log(tau)/2). To see this, split the z integral at1/2.
  Above1/2, z/A<=1/z and its integral is less than1. Below1/2,
  A>=z^2+tau*h/4, giving one-half log(1+1/(tau*h)).
  Since tau*h<=1/4, log(1+tau*h)<=1/4, and
  integral_0^1[-log h]dv=2, the stated bound follows.
- For z^2 and no q, let A0=z^2. Since0<=EP<=1/4 and A>=A0,
  z^2[A0^(-1+EP)-A^(-1+EP)]
  <=(A-A0)/A<=tau/(4A).
  Its integral is at most3sqrt(tau). Higher z powers only improve it.
  The coefficient's own small-q difference contributes a bounded
  multiple of sqrt(tau) because z^2 A0^(-1+EP)=z^(2EP)<=1.

Every rank-zero difference therefore tends to zero uniformly. Finite
sums with bounded coefficients preserve the result; no numerical
assumption about a whole-amplitude order-one prefactor is made.

## Rank-two and rank-four UV-subtracted tensor terms

Here0<=A0<=A<=1. The rank-two Gamma(-EP) term has only a local UV pole.
For its nonlocal difference, the elementary mean-value identity gives

abs((A^EP-A0^EP)/EP)<=log(A/A0).

Since A/A0<=1+tau/(4z^2), its integrated bound is
integral_0^1 log(1+tau/(4z^2))dz<=pi*sqrt(tau)/2.
One proof extends the positive integral to the half-line and
differentiates with respect to b=sqrt(tau)/2:
integral_0^infinity 2b/(z^2+b^2)dz=pi, with zero value at b0.
The exact finite-interval antiderivative is also checked.

For rank four, after removing the local A/EP pole the nonanalytic
function is f_EP(A)=A*(A^EP-1)/EP, up to bounded Gamma/D factors.
Its derivative obeys abs(fprime)<=1+abs(log(A0)) between A0 and A.
Thus its integrated difference is at most
tau/4*integral_0^1[1-2log z]dz=3tau/4.
Finite coefficient differences are controlled by the integrable
A0-log majorants. Gamma factors and their first derivatives are bounded
near EP0. Their finite constants and any retained analytic local finite
counterterms do not obstruct continuity. This is not a prescription
setting those counterterms to zero.

These inequalities extend the UV-subtracted expressions continuously
to EP0 with their explicit logarithms. Accordingly the finite Laurent
coefficient has the same zero-transfer limit, not merely an
unexpanded fixed-EP function.

## The remaining proper graphs

For the one-h two-Phi triangle, group the two massive parameters into
ell and v. The denominator is
ell^2[mu-t(1-v^2)/4]. The raw radial factor is1/(2EP), with the overall
NEGATIVE triangle sign fixed in soft.md. Every tensor numerator power
only shifts the radial exponent by a nonnegative integer; the angular
kernel has positive gap for abs(t)<4mu. The IR-pole and finite
coefficients are analytic on this disk after local UV subtraction.

For the two-h bubble the numerator has two soft derivatives. Tensor
integration gives only q^2 times its scalar nonlocal term, with local
polynomials as usual. Write
Braw=-H(EP)*tau^EP/EP,
H=Gamma(1-EP)(4pi)^(-EP) integral_0^1[x(1-x)]^EP dx, H(0)=1.
Differentiate the positive double-integral representation of H.
The absolute derivative is bounded by products of

integral e^-r max(1,r^-1/4)dr<2,
integral e^-r max(1,r^-1/4)abs(log r)dr<3,
integral_0^1[-log x-log(1-x)+log4pi]dx<5.

For the logarithmic r bound, the pieces are16/9 and at most2/e.
The explicit real-domain primitives and endpoints are checked; no
uncontrolled complex antiderivative branch is used.
Hence abs(Hprime)<13 and
abs(Braw+1/EP)<=13-log(tau).
The vertex contribution tau*(13-log(tau)) tends to zero uniformly.

Both one-h one-Phi seagull bubbles have on-shell denominators depending
only on p^2=mu or pprime^2=mu, not t. Transfer appears polynomially in
their numerators. The h tadpole is scaleless. All contacts, orientations
and the complete EH numerator are thus included.

The only remaining finite-transfer infrared pole is the analytic
one-h triangle coefficient treated in soft.md. Its coefficient and
the finite coefficient are continuous separately. The proof does not
make the finite-transfer IR pole disappear, fix a slope or control an
inclusive scattering experiment.
