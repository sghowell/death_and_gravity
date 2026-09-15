# Entire known formal real/virtual threshold pair

## Conserved current and pole coefficient

At the heavy shell let P=p1+p2, P^2=n and p_i^2=mu. For an outgoing
soft metric momentum q the complete matter current is
J_mn=p1_m*p1_n/(p1.q)+p2_m*p2_n/(p2.q)-P_m*P_n/(P.q).
It obeys q.J=0 by momentum conservation. In the heavy rest frame its
harmonic contraction equals its physical transverse-traceless sum.
Keeping every self and cross term, the angular average of
omega^2[J:J-(trace J)^2/2] is
R3=mu+n/2+2V*I0/n-2(n-2mu)*I0
  =mu+n/2-4mu(1-mu/n)*I0,
I0=atanh(beta)/beta, V=n^2-4mu*n+2mu^2.
The two light self terms give mu, the heavy self term n/2,
the light-light cross term2V*I0/n, and both heavy-light cross terms
-2(n-2mu)*I0. None is dropped before forming the conserved current.

The pointwise kernel is Q^2/(2n)*(1-x^2)^2/(1-beta^2*x^2)^2,
Q=n-4mu. Let F0 be its normalized angular integral without Q^2/(2n).
Then R3=Q^2*F0/(2n)>0. This matches the soft limit of the complete
four-graph S6.291 production tree, not a single emission graph.

The three virtual pair triangles contribute real pole
+4mu(1-mu/n)*I0/(16pi^2 kappa e).
All external sqrtZ factors add-(mu+n/2)/(16pi^2 kappa e).
Hence Re(deltaT/g)=-R3/(16pi^2 kappa e)+finite.
The imaginary+i*pi*V/(16pi^2 kappa*n*beta*e) remains.
Using S6.291 K0(n)=g^2*R3/(4pi*kappa),
K0(n)/(2e)+pi*g^2*2Re(deltaT/g)_pole=0 exactly.

This algebra concerns the real heavy delta coefficient only. The
imaginary Coulomb contribution multiplies the heavy principal value
and other cuts. It has not been canceled or declared irrelevant.

## Whole finite compact-window functional

Keep a fixed C1 function W on[n,n(1+L)],0<L<=1/2. S6.291 supplies
rho_e(s)=(s-n)^(-1+2e)*K_e(s).
Use y=(s-n)/n and H_e(y)=n^(2e)*K_e(n(1+y))*W(n(1+y)).
Its exact subtraction-and-addition formula is
integral_0^L y^(-1+2e)*H_e(y)dy
=H_e(0)*L^(2e)/(2e)
 +integral_0^L y^(-1+2e)[H_e(y)-H_e(0)]dy.
Differentiability and the massive angular gap justify the finite
expansion, while W being fixed rules out hidden e-dependent weights.

With F0=integral_0^1(1-x^2)^2/(1-beta^2*x^2)^2 dx and
L_B=integral_0^1(1-x^2)^2*ln(1-x^2)/(1-beta^2*x^2)^2 dx,
the complete S6.291 first coefficient at threshold is
H1(0)=K0(n)W(n)*
 [EulerGamma+1+ln(n/(16pi nu^2))+L_B/F0].
Its finite local contribution is H1(0)/2+H0(0)ln L.

The complete known virtual finite coefficient contributes
K0(n)W(n)/2*
 [Re(B1)/R3-EulerGamma+ln(4pi nu^2)].
Combining these with the exact regular integral yields
lim_(e->0+)D_W(e)=
 integral_0^L[K0(n(1+y))W(n(1+y))-K0(n)W(n)]/y dy
 +K0(n)W(n)ln L
 +K0(n)W(n)/2*[1+ln(n/4)+L_B/F0+Re(B1)/R3].
This uses B1 from the whole D tensor/master sum, including all
evanescent terms. Neither a finite delta counterterm nor a numerical
matching constant has been selected to obtain it.

EulerGamma and nu^2 cancel exactly. For a common rescaling of all mass
squares by lambda, B_e rescales as lambda^(1+e). Since Re(B0)=-R3,
Re(B1)/R3 gains-ln lambda, canceling the logarithm of n.
Equivalently set Bhat_e=B_e/n^e. Then
Re(Bhat1)/R3=Re(B1)/R3+ln n, and the local bracket reads
1-ln4+L_B/F0+Re(Bhat1)/R3. The whole coefficient is unit consistent.

## Explicit conditional remainder

In the notation of S6.291 define on0<=e<=1/8,0<=y<=L
M2=sup|d_e^2[H_e(0)*L^(2e)]|,
M10=sup|d_y H_e(y)|,
M11=sup|d_e d_y H_e(y)|.
Its complete real remainder is bounded by
e*[M2/4+2M10*L*(1-ln L)+M11*L].

For the virtual part put
C_e=Gamma(1-e)*(4pi nu^2)^(-e)*B_e,
Mvirt2=sup_(0<=e<=1/8)|Re d_e^2 C_e|.
The physical-root subtraction in notes/masters.md, positive mass
gaps and regular Gamma factor prove this supremum finite for each
specified positive massive domain. Taylor's theorem applied to C_e/e
adds at most
e*g^2*abs(W(n))*Mvirt2/(16pi*kappa).
The sum is an explicit conditional O(e) bound for the named whole
functional. It is not an evaluated detector error or a uniform bound
as beta approaches an excluded massless/threshold endpoint.

Independent whole-D tests evaluate the real cut and virtual expression
before subtraction, then compare with this finite prediction at
decreasing positive e, for two unequal-mass/scale cases. These tests
check the signs, normalization, scale and O(e) behavior. They are
finite numerical checks, not formalized all-domain theorems.

## Physical interpretation boundary

D_W is a deliberately named partial distributional combination, not the
complete amplitude or a positive physical spectral measure. The
unstable heavy resonance requires additional uniform/width analysis.
The H-metric response and local/higher-EFT matching contribute finite
pieces outside this minimal proper-plus-LSZ coefficient. Coulomb
principal values, other cuts, a fully specified physical infrared
observable, fixed-transfer Regge data and original P8 remain open.
