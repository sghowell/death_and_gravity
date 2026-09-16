# Complete inherited finite graph assembly and matching separation

With channels(s,t,u), q_a=n-a, V(a)=(a-2)^2-2 and
Htree=sum1/q_a, use c_E=EulerGamma-ln(4pi E^2). The known S294
nonendpoint bracket is

F_g=Htree{sum[V J1+2J0]+c_E[sum V J0-2]}
 -sum K1(a)+sum[4(a-2)T(a)-4(a^2-2n)U(a)/q_a]
 +2sum_(a!=b)V(b)W(a,b)/q_a,

where K1,T,U are the whole frozen masters and W is the exact collapsed
box in notes/masters.md. The code replaces each of the six original K
integrals symbolically and checks equality with the entire S294 bracket.

For the full S293 quartic sector the raw expression is
C Gamma(-e)(4pi nu^2)^(-e)[-2T_e+E_e]/(16pi^2 kappa).
At original mu=nu^2=1 its complete first coefficients are

T0=sum VJ0/2-1,
T1=sum[VJ1/2+J0+(a-2)Lll]-1,
E0=5/3,
E1=-sum(a+2)Bbar-1.

The last -1 in E1 is the full evanescent constant; it is not discarded
because it has zero forward second derivative. The local UV pole is
-E0/e and the gravitational IR pole is2T0/e in these stripped units.
Remove the local UV pole first, then subtract exactly the inherited
2T0(E^2)^e/e soft reference. With c0=EulerGamma-ln(4pi), the named
pure-pole known finite representative is

F_C=2T1-E1+2c_E T0-c0 E0.

The code checks this Laurent expansion before any estimates. A finite
local counterterm remains separate: this formula does not select its
physical value. The auxiliary E is not a detector-resolution choice.

## Whole endpoints and no double counting

For N_a=2-2a-bc, the complete harmonic contraction of both corrected
vertices produces

A_endpoint=-sum[(2N_a+a+a^2/2)f1(a)/a+(a+2)f2_g(a)]/kappa.

Use the whole S290 f1 with its OS subtraction and f2_generated with
quartic=0. Both massive triangle assignments, the entire kinetic
term and transverse H response remain. Setting only this endpoint's
C argument to0 prevents adding S293's direct C bubble twice; it does
not set the original contact to0. The complete selected known sum is

A_known=[g^2 F_g+C F_C]/(16pi^2 kappa)+A_endpoint.

The constant curvature and RH coordinates would contribute
-10ell/kappa-h sum(a+2)/(kappa(n-a)); an independent finite local
constant is also left explicit. No value is assigned. These displayed
minimal axes do not exhaust higher-derivative or other-sector matching.
Changing them can change the matched physical rate, so the known bound
must never be labeled a bound on arbitrary matching.

The matter-only loop, pure Einstein/Phi loop and spectator metric-loop
sectors have their own owners and are not included again here. Real
radiation, the physical soft conversion and all-resolution leading-soft
sum are also distinct from the bounded hard representative.
