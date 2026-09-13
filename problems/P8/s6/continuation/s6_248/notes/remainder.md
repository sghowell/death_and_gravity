# Strong normal form for the actual two-species reference

Fix a finite external radius Pmax and a finite number j of simultaneous time derivatives. Both physical masses are fixed. Constants below may depend on them, the whole smooth reference and Pmax; neither uniform mass nor analytic-in-j estimates are claimed.

## Same actual scalar state with arbitrarily deep comparisons

For the minimally coupled scalar put chi=a^(3/2)phi. Its exact mode equation is chi''+omega_eff² chi=0, with omega_eff²=k²/a²+n-(3/2)H'-(9/4)H². On a compact time interval, sufficiently large k makes every specified finite WKB iteration positive. The iteration
W_next²=omega_eff²-(1/2)W''/W+(3/4)(W'/W)²
has a full finite inverse-k symbol expansion with smooth coefficient jets. Each deeper iteration reduces the residual by further inverse powers. Its exact KG solution initialized with the full WKB Cauchy pair at-1 differs by variation of constants by an arbitrarily prescribed inverse power, after allowing for the finite number of time derivatives. The frequency energy estimate is uniform on that compact high-k region.

The actual S240 SLE is defined by the original compact smooth sampling bump and is basis independent. Express that SAME state in a sufficiently deep exact comparison basis. The diagonal sampled energy is bounded above and below by positive constants times k. Every off-diagonal WKB term has nonstationary phase derivative2W comparable to k. Repeated integration by parts through the full smooth bump has no boundary term and supplies any prescribed inverse power; the exact-comparison error can be made smaller by deepening the iteration first. The SLE Bogoliubov formula and its strict energy cone then give the same arbitrarily high algebraic smallness for mixing. Exact mode time derivatives cost only finitely many k powers, controlled by selecting the depth AFTER j and the desired decay have been fixed.

This is also the applicable scalar result behind [Olbermann, arXiv:0704.2986v2, Lemmas4.4-4.5 and Theorem4.9](https://arxiv.org/abs/0704.2986v2): deep adiabatic comparisons approach the SLE to arbitrary fixed regularity, and the SLE is Hadamard. The paper does not give the coupled inverse, the constants here or a finite-dimensional continuation of the actual state. Those require the preceding mode argument and the model-specific assembly below.

The original Proca all-order prepared mode comparison used in S228/S229 remains part of the unchanged reference. The new scalar argument extends the actual second species; it is not a claim that the old QG1 inverse already applies to the new summed operator.

## Entire singular and finite matching before the remainder

For ellvec=Pvec-kvec, split the ORIGINAL internal integral into a bounded region and k>R, where R>2Pmax and exceeds the fixed mass/background thresholds. This is an estimate, not a physical cutoff. Both legs are large in the second region, uniformly over the fixed external ball.

The complete high-radius scalar stress amplitude is (trQ-nhat Q nhat)/2. Its full dimension invariants are
C_H(d)=1/[2d(d+2)], B_H(d)=(d²-3)/[4d(d+2)].
Add the full Proca invariants
C_P(d)=(2d²+d-8)/[2d(d+2)],
B_P(d)=(d³-4d²-3d+16)/[4d(d+2)]
BEFORE fixed-physical contraction with Q=2wI_3-2cPi.

The resulting matrix M(d)=(C_H+C_P)[[12,-4],[-4,4]]+(B_H+B_P)[[36,-12],[-12,4]] has
M(3)=[[8,-8/3],[-8/3,32/15]]
 =L0^T diag(8,56/45)L0,
M'(3)=[[68/15,-68/45],[-68/45,224/225]],
L0=[[1,-1/3],[0,-1]].
S247 independently derives these from the full stress cuts and finite actions. The full same-prescription fourth local matrix is
-L0^T diag(2(log(n)+2),(2/45)(log(n)+2))L0.
No massless/log-only substitute, altered subtraction scale or omitted heavy cut is used.

At frozen scale a0, a0^(-d-2)k^(d+1)dk becomes p^(d+1)dp under k=a0p for EVERY nearby d. This all-dimensional identity and the physical first jet fix the highest finite contact, not merely its pole. The undifferentiated dimensional finite part is the original S245 whole W6 continuation with absolutely convergent actual-state subtraction. Take that fixed limit first; deepen the physicald3 comparison only to establish higher time regularity. Do not invent an all-dimensional SLE by analytic continuation.

With normalized output64pi²/a(t)³, the full leading radial matrix is
32 M(3) k^4 sin(2k Delta_sigma)/[a(t)^4 a(s)].
Its Abel integral is24M(3)/[a(t)^4 a(s)Delta_sigma^5].
The proper-time ratio tau^5/[a(t)^4 a(s)Delta_sigma^5] is smooth at tau0 and equals1 there. Its first lag coefficient is-3a'(s)/(2a(s)), explicitly checked.

At finite transfer,
omega_k+omega_ell=2k/a-(nhat dot P)/a+O(k^-1).
KEEP exp(-i(nhat dot P)Delta_sigma) as an orderzero high-k amplitude. It equals1 on the diagonal and its difference contains a lag. Misclassifying it as inverse-k-small is false. Both it and the proper-time ratio therefore reduce the unmatched leading singular order by at least one.

For each simultaneous derivative Dt+Ds, the leading phase contributes k[a(t)^-1-a(s)^-1], which has a compensating lag. Repeated derivatives preserve conormal order: every extra radial power comes with enough lag factors. The finite-transfer phase and smooth geometric amplitudes have bounded fixed derivatives uniformly on the ball. A sufficiently deep two-leg comparison leaves an absolutely integrable error even after the specified jets. On bounded internal momenta, the full massive covariance and smooth time evolution give ordinary smooth remainders, including zero internal legs. Complete physical sums, not separate singular polarization charts, are bounded there.

## Four OUTPUT primitives

After the matched complete q0 reference is removed, all remaining finite-part singularities have order at most4. Highest fourth local terms have already been matched; remaining local terms are at most third order in the spatial quantum block. The two full Ward legs, distinct clock contacts, fixed profiles and classical terms are local of the mixed orders established in assembly.md.

A causal finite part tau^-n, n<=4, has fourth primitive proportional to tau^(4-n)log(tau), plus the polynomial fixed by its original extension. Expand a smooth two-time amplitude about the diagonal to the finite depth needed to make the remainder ordinary. Simultaneous derivatives obey the same estimate. All original finite extensions and upper contacts stay in this operation; preparation only removes the unchanged LOWER zero germ.

For a local term c_j(t)D^j, its n-fold OUTPUT primitive has integral kernel
sum_r (-1)^r binom(j,r) (t-s)^(n-1-j+r)c_j^(r)(s)/(n-1-j+r)!.
For j<n the sum starts at0. For j=n it starts at1 and has the additional multiplication c_n(t). Thus the eta pivot gives -6delta(t)² multiplication PLUS its complete lower kernel, not a freely commuted constant. Code checks every local case for row orders2 and4.

The resulting ACTUAL full normal form is
I_rows T_adapt=Fref+V,
I_rows=diag(I4,I4,I4,I2),
Fref=diag(-6delta(t)²,gamma L0^T diag(Ftrace,total,(8/3)F2,total)L0,-1),
gamma=1/(64pi²kappa).

For every fixed j,
ess sup_(0<|P|<=Pmax) max_i sum_l |(Dt+Ds)^j V_il(t,s;P)|
 <=C_j(Pmax)(1+|log(t-s)|).
Each C_j is finite and contains all current coefficients, masses, finite contacts, actual state terms and normalization. They are NOT evaluated numerically. This strong statement follows from the full actual mode and singular-extension analysis, not the derivative-losing weak response estimate.
