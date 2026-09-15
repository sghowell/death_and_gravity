# Literal two massive insertions and complete generated F1

For the triangle with active scalar mass squared a and spectator b,
write the internal spectator momentum k and the two active momenta
r=p-k and rprime=pprime-k. The exact tensor identity is
q.T=(rprime^2-a)r-(r^2-a)rprime.
All scalar momenta are kept off shell when using this identity.

The normalized triangle master has the opposite sign to the bubble
master. Summing both line insertions and translating the common mixed
bubble gives pprime B(pprime^2)-p B(p^2). Each of the two literal metric
cubic contacts contributes-eta B, so their combined contraction is
-q[B(pprime^2)+B(p^2)]. The total is exactly
p B(pprime^2)-pprime B(p^2), the mixed self-energy Ward identity.
Signs are independently checked from the ordinary+ig, metric-iT and
metric-cubic+ig eta factors and all propagator/master factors.

These massive identities may first be evaluated in a Euclidean
absolutely convergent dimensional strip, with the complete tensor
integrals and metric contacts present, and then meromorphically
continued. No massless infrared assumption enters these graphs.

On the external mass shells p^2=pprime^2=mu, set
x=(1-z)(1+v)/2,y=(1-z)(1-v)/2,0<=z,v<=1.
The combined simplex measure is(1-z)dz dv and
Delta_ab=(1-z)a+zb-z(1-z)mu-xy*t.

Translate k by xp+ypprime. The PP part of the literal stress tensor is
2z^2 PP. An independent02 projection in a spacelike-transfer frame has
no eta or qq component; odd loop terms and the off-diagonal isotropic
ell0*elly moment vanish. Thus, with c=g^2/(16pi^2), the whole generated
renormalized correction is

f1(t)=c sum_(a,b)=(mu,n),(n,mu)
 int_0^1dz int_0^1dv (1-z)z^2[1/Delta_ab(t)-1/Delta_ab(0)].

The subtraction is the inherited kinetic counterterm. Reversing z in
the heavy-line t0 integral makes the two integrands sum to
z(1-z)/F(z), F=mu(1-z)^2+nz: exactly the complete mixed Pi_prime
integrand. Therefore F1(0)=1 at this formal loop order.

Differentiation before integration is controlled by the full positive
denominators. The v integration and the same z reversal give
f1prime(0)=c/6 int z^2(1-z)^2/F^2 dz=Pi_second(mu)/6>0.
At every higher order k>=1 the coefficient is
int_0^1 (1-v^2)^k dv/4^k times
sum_(a,b) z^2(1-z)^(2k+1)/Delta_ab(0)^(k+1).
An independent simplex beta integral checks the normalization.

Metric contact, quartic bubble/tadpole and H-metric mixing tensors have
no PP component, so do not alter this F1 projection. They remain
necessary for a full F2/physical endpoint calculation. These formulas
also do not forbid an additional finite Ricci-derivative local term.
