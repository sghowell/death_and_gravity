# Complete production tensor and full physical polarization sum

Take p1^2=p2^2=mu,q^2=0,P=p1+p2,R=P-q,R^2=n,s=P^2 and
T_m(p,r)=p_mu r_nu+r_mu p_nu-eta_mu_nu(p.r-m).
The original propagator/vertex products give, after stripping g/sqrt(kappa),

A=T_mu(p1,p1-q)/D1+T_mu(p2,p2-q)/D2+T_n(P,R)/(s-n)+eta,
D_i=(p_i-q)^2-mu.

The literal phase product (ig)(-iT)i=+igT and cubic contact+ig eta
fix the displayed relative signs. The contractions withq are respectively
-p1,-p2,+R,+q. They sum0. Deleting the contact gives-q; deleting heavy
emission gives-R. Their TT projections vanish, but their Ward role
cannot be discarded.

In COM choose p1=(E,a,0,b),p2=(E,-a,0,-b),q=(omega,0,0,omega),
R=(2E-omega,0,0,-omega),mu=E^2-a^2-b^2,n=4E(E-omega).
D1=-2omega(E-b),D2=-2omega(E+b),s-n=4Eomega.
For epsilon=(0,1,+/-i,0)/sqrt2, the complete single-helicity amplitude is
-g E a^2/[sqrt(kappa)omega(E^2-b^2)].
Writing Q=s-4mu,B=1-4mu/s,x=cos(theta), it is

A_hel=-g Q(1-x^2)/[sqrt(kappa)(s-n)(1-Bx^2)].

The general-D sew follows from the FULL physical projector in d=D-2
transverse dimensions:
P_ij,kl=(delta_ik delta_jl+delta_il delta_jk)/2-delta_ij delta_kl/d.
Delta contraction proves idempotence, tracelessness and
rank=d(d+1)/2-1=D(D-3)/2. The isotropic transverse part of A is killed;
its remaining rank-one tensor contracts to(a^2)^2(1-1/d). Thus the whole
physical sum is4(D-3)/(D-2)*g^2 Q^2/(kappa(s-n)^2)*
(1-x^2)^2/(1-Bx^2)^2. Literal component matrices in D4/5/6 and independent
rational COM fixtures check every projector and Ward entry.
The general D proof is the displayed delta contraction, not extrapolation
from those three integer dimensions.

At D4,H and h are DISTINCT. Optical half of phase beta_Hh/(8pi), not
the identical-pair quarter, gives
rho_0=g^2 Q^2 F0(B)/[8pi kappa s(s-n)],
F0(B)=integral0^1(1-x^2)^2/(1-Bx^2)^2 dx.
The full rational decomposition is
1/B^2+2(B-1)/[B^2(1-Bx^2)]+(B-1)^2/[B^2(1-Bx^2)^2].
With I0=atanh(sqrtB)/sqrtB and I2=1/[2(1-B)]+I0/2 this yields
F0=[3-B+(B-1)(B+3)I0]/(2B^2).
Both endpoints are removable: F0(0)=8/15,F0(1)=1.
Pointwise bounds and positive B derivative prove8/15<=F0<=1.

This is the complete production amplitude and its FORWARD sew. It does
not establish the nonforward mixed cut, exact heavy stability or finite G.
