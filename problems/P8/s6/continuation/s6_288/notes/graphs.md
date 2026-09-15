# Complete minimal graph inventory and physical pole theorem

## Finite topology count

Let M be the number of minimal two-Phi vertices and G the number of
pure-graviton vertices in a connected four-Phi one-loop graph without
external gravitons. There are Iphi=M-2 scalar lines and Ih=G+2
graviton lines. Each matter vertex has at least one graviton end, and
each pure-gravity vertex at least three. Consequently M+G<=4, M>=2.
Writing e for the total pure-gravity valence above3G gives exactly

 (M,G,e,matter_h_ends,Iphi,Ih)=
 (2,0,0,4,0,2), (2,1,0,3,0,3), (2,1,1,2,0,3),
 (2,2,0,2,0,4), (3,0,0,4,1,2), (3,1,0,3,1,3),
 (4,0,0,4,2,2).

These count all internal valences, not a selected helicity set.
For G0,M4 the four cubic matter vertices give the alternating boxes,
proper one-h/two-Phi vertices with an exchanged h, Gaussian Phi
metric bubbles, or metric one-point attachments. For G0,M3 the
h-valences2,1,1 give genuine two-h/one-Phi seagull triangles, proper
seagull bubbles with an exchanged h, or Gaussian metric tadpoles.
For G0,M2 the valences2,2 give genuine two-h bubbles, and3,1 give
a scaleless h tadpole at a matter contact.

The G1 cubic classes give the entire EH proper vertex triangle or
EH/seagull bubble, or one-point attachments. G1 with a quartic
metric vertex gives the scaleless pure-h metric tadpole. G2 gives
the massless metric bubble. Ghost loops are metric-reducible in
this one-loop four-scalar count: one attachment is a one-point
function, two supply a metric self-energy, and extra independent
attachments to the two scalar paths add a loop. They are included
in that whole metric kernel. The scalar external-leg loops are
included separately through the complete S283 residue derivative.
Pure scalar higher-arity/source vertices are not silently counted
as minimal vertices; their distinct sectors are not claimed absent.

## Literal full tensor-box reduction

Set p3=p+q,p4=r-q, p^2=r^2=mu, q^2=t,
p.q=-t/2,r.q=t/2 and b=(p+r)^2. Route the denominators as
d1=k^2, d2=(p-k)^2-mu, d3=(k+q)^2, d4=(r+k)^2-mu.
The general-D harmonic contraction of two scalar stresses is

 2[(x.y)(v.w)+(x.w)(v.y)]-2(x.v)(y.w)
 +2mu[(x.v)+(y.w)]-2Dmu^2/(D-2).

Direct substitution at both exchanged gravitons yields the ENTIRE
box numerator, with z=b-2mu:

 [VD(b)+z(d1-d2-d4)] [VD(b)+z(d3-d2-d4)].

Independent full symmetric-tensor contractions in D4,5,6 confirm
this identity with all loop components arbitrary. The difference
from VD(b)^2 is

 VD*z*(d1+d3-2d2-2d4)
 +z^2*[d1d3-(d1+d3)(d2+d4)+d2^2+2d2d4+d4^2].

Every term cancels at least one propagator. The remainder is
C0mumu(b), C00mu(t), Bmm(b), on-shell mixed bubbles, B00(t), and
at most rank-two C00mu tensor moments. No external1/t factor or
Gram division was introduced. Thus subtracting the known scalar
box is now justified for the pole comparison; it was NOT assumed
identical to the whole tensor diagram.

## Transfer-weighted remainder bounds

First choose a compact subthreshold spectator interval
delta<=b<=4mu-delta and0<tau=-t<=delta/2, with a smaller tau cap
if needed. Then c=4mu+tau-b also lies below4mu, and both massive
box/triangle angular gaps are at least delta/8.
A bounded complex Breit family realizes these on-shell invariants:
E=sqrt(b)/2, transverse component
i*sqrt(mu+tau/4-b/4), longitudinal components+/-sqrt(tau)/2.
Its coefficients are smooth in sqrt(tau) and uniformly bounded.
Parameter denominators of the transfer triangles are positive;
no physical-threshold boundary limit is interchanged in this step.

S287 gives, after rescaling the mass unit,
tau*abs(C00mu_raw)<=C*sqrt(tau) for pole/finite coefficients.
For rank-one or rank-two moments shift the loop before integrating.
The resulting finite numerator has bounded coefficients and the
same positive scalar denominator. Odd shifted moments vanish.
Rank-two UV-subtracted moments contain integrable log denominators,
bounded by a constant times1+abs(log z) on the scalar parameter
interval; they have no1/t singularity. All regular UV pole
coefficients are polynomial and also vanish after multiplication
by t. The fixed-spectator C0mumu(b), Bmm(b) and on-shell Bmix
pole/finite coefficients are bounded in tau. B00(t) grows only
as log(tau). Therefore t times the ENTIRE denominator-cancelled
box remainder tends to zero. The same argument covers the genuine
seagull triangle and bubble; their full numerators have ranks at
most2 and polynomial external coefficients.

## Whole proper stress tensor, including possible singular F2

Using only the F1 projection would be insufficient. Apply S287's
same numerator-ideal and parameter estimates to every unprojected
component instead. Its EH triangle has rank<=4 and soft ideal
(k,q)^2; after the shift the L0 part is in(z,q)^2. The established
q^2, zq, z^2 and UV-subtracted L2/L4 bounds are independent of
which external tensor component multiplies a monomial. All of
those finite external coefficients are bounded in the family
above. The other proper topologies have the same gapped analytic
or q^2*log(tau) bounds already proved there.

Thus the complete proper stress tensor, not only F1, is continuous
in its pole and finite coefficients. The zero-q ordinary Ward
bridge and Sigma(mu)=0 give Gamma_loop(0)=2pp*Sigma_prime.
The two endpoint scalar factors give-Sigma_prime*Gamma_tree;
their entire constant tensor cancels. Consequently t times the
complete vertex-corrected exchanged amplitude tends to zero.

In Gamma=2PP F1+(qq-eta*t)F2 this statement controls the ENTIRE
t*F2 tensor. It does not assume F2 finite, discard it, or prove a
finite F1 slope. A possible1/sqrt(-t) form factor is compatible
with the proved tensor continuity and cannot generate a simple
exchange pole in this argument.

## Metric reducible pieces and analytic continuation

The massless graviton/ghost metric two-point kernel has four soft
derivatives and no mass scale, hence its nonlocal term is
t^2*log(-t), modulo local degree-four terms. Scaleless metric
tadpoles vanish in the retained dimensional scheme. Between two
conserved scalar currents its two propagators cancel those t^2
factors; t times the result vanishes.

The full gapped Phi Gaussian metric kernel includes both bubble
and metric-contact terms. S285's matched whole volume density
sets its constant metric term to zero, removing the double pole.
The only residual simple pole is the Newton shift
delta_kappa_Phi=2c_Phi_R. Quadratic-curvature and higher-subtracted
gapped terms are regular in this pole comparison. The full
Newton-shape sign is obtained by expanding
A_tree(kappa+delta_kappa), not guessed from a curvature convention.

We have therefore determined the rational principal parts on
an open subthreshold spectator interval after retaining all scalar
boxes explicitly. A rational identity holding there continues
to generic complex invariants, with physical boundary values
taken afterward. This does not exchange a physical forward and
IR limit, promise uniformity at the massive threshold, or convert
the Coulomb branch into an isolated meromorphic pole.
