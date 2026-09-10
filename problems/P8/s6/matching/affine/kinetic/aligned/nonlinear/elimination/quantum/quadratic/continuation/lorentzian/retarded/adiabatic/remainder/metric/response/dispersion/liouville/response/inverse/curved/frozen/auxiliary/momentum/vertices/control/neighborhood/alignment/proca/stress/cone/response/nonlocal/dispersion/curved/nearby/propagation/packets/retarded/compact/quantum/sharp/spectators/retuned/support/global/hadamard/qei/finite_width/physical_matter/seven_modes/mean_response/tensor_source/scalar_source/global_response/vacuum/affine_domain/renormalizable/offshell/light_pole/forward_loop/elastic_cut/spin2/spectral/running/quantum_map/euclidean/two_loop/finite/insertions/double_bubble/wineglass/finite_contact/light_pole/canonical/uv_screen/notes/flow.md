# Explicit new candidate and one-loop marginal flow

Use canonical gauge kinetic term -F^a F^a/4, a=g_s^2 and
Q=16 pi^2. SU(3) generators have Tr(Ta Tb)=delta_ab/2
and Casimir 4/3. The code constructs all eight matrices and
checks every trace, Hermiticity and the Casimir sum.

There are fourteen fundamental four-component Dirac flavors.
They are vectorlike, so left/right gauge anomalies cancel.
Only the pair q+ and q- has Yukawa interaction
-y Phi (bar(q+)q+ - bar(q-)q-). Their common Dirac mass mF
preserves Phi parity by exchanging q+ with q-. This is not
a chiral symmetry which the mass term would break. The twelve
other Dirac flavors have no Phi Yukawa. The scalar H remains
a singlet without a direct fermion coupling.

The one-loop renormalization formulas identified in the literature
note are inputs, not newly proved loop integrals. For the six
active color/flavor entries set Y=diag(y,y,y,-y,-y,-y).
The code evaluates their matrix products and traces before
reducing to one Yukawa. It independently differentiates V=L Phi^4/24
in the scalar potential tensor formula. With t=log(mu/muF):

    Q da/dt = -(10/3) a^2,
    Q dy/dt = y (15 y^2 - 8a),
    Q dL/dt = 3L^2 + 48y^2 L - 288y^4.

The neutral, dimension-one H Phi^2 interaction cannot enter
these one-loop marginal logarithms: its extra inverse momentum
power makes such marginal contributions superficially convergent.
This statement does not remove relevant-parameter counterterms.

Set y^2=ry a, L=rL a. The Yukawa equation gives ry=19/45.
The quartic equation is

    675 rL^2 + 5310 rL - 11552 = 0.

Its derivative is positive for rL>=0. Its values at 1 and 2
have opposite signs, so there is exactly one positive root
between them: rL=(sqrt(65985)-177)/45. The other root is
negative and is not an admissible positive quartic ray.

For arbitrary a0>0 and t>=0, a=a0/[1+(10/3)a0 t/Q]
is positive and decreases to zero. Together with the two
ratios it solves all three displayed equations; the code
differentiates each solution and checks its initial value.

This is exact one-loop marginal flow, not a full UV verdict.
Relevant couplings, fermion thresholds, the quantum potential,
nonperturbative dynamics and bounce matching remain separate.
