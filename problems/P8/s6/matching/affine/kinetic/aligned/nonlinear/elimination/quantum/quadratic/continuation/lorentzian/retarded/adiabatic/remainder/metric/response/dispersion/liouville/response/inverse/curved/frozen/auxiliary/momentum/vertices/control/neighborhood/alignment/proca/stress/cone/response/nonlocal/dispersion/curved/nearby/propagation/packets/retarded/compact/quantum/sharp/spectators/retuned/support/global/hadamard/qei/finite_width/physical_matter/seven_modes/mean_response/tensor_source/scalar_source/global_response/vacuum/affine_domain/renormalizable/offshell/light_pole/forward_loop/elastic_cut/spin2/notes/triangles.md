# Actual mixed triangles and the Ward-normalized slope

Write M=mu^2 for the heavy mass squared and g=G^2. For a triangle
whose two stress-adjacent lines have mass squared mA and other line
has mass squared mB, choose Euclidean routes k, k+q and k-p.
Let chi be the single-line parameter and z the second equal-line
parameter. The literal four-dimensional shift is l=k+zq-chi p.

Its parameter denominator before the on-shell restriction is

Delta=(1-chi)mA+chi mB+z(1-z)q_E^2
      +chi(1-chi)p_E^2+2z chi p_E.q_E.

With p_E^2=-1, q_E^2=-t and p_E.q_E=t/2, this becomes
(1-chi)mA+chi mB-chi(1-chi)-z(1-chi-z)t.
The full vector shift and both mass assignments are checked exactly.

Under the null projection e.q=0, e.k=e.l+chi e.p.
Odd loop terms vanish and the tensor trace is proportional to e^2=0.
The surviving stress numerator is 2chi^2(e.P)^2. The three-denominator
Feynman factor is two; the finite radial integral is
integral d^4l/(2pi)^4 (l^2+Delta)^(-3)=1/(32pi^2 Delta).
Dividing by the external tensor 2(e.P)^2 gives prefactor g/(16pi^2).

Put d1=xM+(1-x)^2. Rescaling the remaining pair parameter to z in [0,1]
gives the two unit-square terms

F_loop(t)=g/(16pi^2) integral_0^1 dx integral_0^1 dz
 [x^2(1-x)/(d1-(1-x)^2 z(1-z)t)
  +x(1-x)^2/(d1-x^2 z(1-z)t)].

The first term inserts stress on a light line, the second on a heavy
line. Both are needed. At zero transfer their weights sum to x(1-x).
Therefore F_loop(0)=Pi'(1), the actual light-residue quantity calculated
independently in S6.112. The same delta Z=-Pi'(1) gives
F_R(t)=1+F_loop(t)-F_loop(0), with F_R(0)=1.

Differentiating and integrating the pair parameter yields

F_R'(0)=g/(96pi^2) integral_0^1 [x(1-x)/d1]^2 dx
       =Pi''(1)/6 >0.

This relation is checked against the actual parent two-point kernel
and couplings, not just a symbolic Ward ansatz.
Using x(1-x)/d1<=(1-x)/M and pi>3,

0<F_R'(0)<g/(288pi^2 M^2)<g/(2592M^2)<10^-405

at the actual rational calibration.

The two normalized transfer weights obey
rho_L=(1-x)^2 z(1-z)/d1<=1/4 and
rho_H=x^2 z(1-z)/d1<=1/(4M)<=1/4.
Positive-gap identities prove these inequalities. Thus F_R is analytic
on |t|<4. For |t|<=r<4 the exact resolvent identity gives

|F_R(t)-1| <=F_R'(0)|t|/(1-r/4),
|F_R(t)-1-tF_R'(0)| <=F_R'(0)|t|^2/[4(1-r/4)].

The majorants are integrable over the parameter square and justify
differentiation and the limiting estimates. The transfer light-pair
threshold is not confused with the mixed two-point threshold.
