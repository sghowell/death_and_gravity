# Direct complex-time jets on the actual nearby solution

Use S6.89's exact QQ rational field in (u,N,z), with z=e(p-b).
The actual S6.88 holomorphic solution on |u|<=T=10^-7 lies in the
whole box with radii (T,4*10^-7,3*10^-5), centered at (0,1+10^-6,0).
S6.90 establishes the required nonzero complex annuli. These are
bounds on the actual solution, not a sample grid or a new ODE.

For a shifted polynomial P=P0+sum_{alpha!=0}P_alpha x^alpha,
let dP=sum |P_alpha| r^alpha. On the full complex polydisc,
|P|<=|P0|+dP and |P|>=|P0|-dP when the latter is positive.
The same coefficient calculation bounds each partial derivative.
For a nonzero rational f=P/Q,
|partial_j f/f| <= sup|partial_j P|/inf|P|
                 +sup|partial_j Q|/inf|Q|.
The implementation checks both nonzero denominators, including
the numerator P in this logarithmic expression. It does not use
ordering of complex numbers. Its rational-box coefficients feed
a modulus annulus explicitly.

The actual on-constraint flow is N'=-FORCE/PIV and
z'=Z0+Kz N'. Complete modulus bounds give
|N'|<7/50000, |z'|<9, |H|<1/62500 and
|e'/e|<1/12500, with |PIV|>2 retained.
The e derivative uses e'/e=Lt+Ln N', not a frozen-lapse derivative.
The code bounds all Cartesian partial logarithmic derivatives of
D, h=(1+u^2)^3, M, PIV and C=c_clock^2. All 15 declared caps
are checked against native whole-box rational enclosures.

Contracting those caps with (1,7/50000,9) bounds total time jets.
In particular |C'/C|<=363/100000.
The literal action root A_c=sqrt(K_c omega_c) satisfies
A_c^2=-D PIV sqrt(C)/(4 N h R M^2), using e^4=N^2 h/D.
Consequently
A_c'/A_c=(D'/D+PIV'/PIV+C'/(2C)-N'/N-h'/h-H-2M'/M)/2.
Its declared upper bound is 1/400. The N term uses |N|>99/100,
also for complex time; real N>1 is not substituted into that step.

For A_m=sqrt(K_m omega_m)=e/sqrt(R),
|A_m'/A_m|<1/10000. The original conserved matter charge is
ell R^3, hence ell'=-3H ell and |ell'|<3/500000.
Finally omega_m=N/(eR), omega_c=sqrt(C) omega_m.
Using both logarithmic flow terms and |omega_j|<4 yields
|omega_j'|<1/100 for both signed pairs.

All roots are the branches positive at the real central point.
Their factors are holomorphic and nonzero on the simply connected
actual time disc, so these branches and logarithmic derivatives
are legitimate there. Fixed-phase derivatives and the old
primitive terms have already been retained in the pinned
rational flow; this calculation differentiates that actual flow.
