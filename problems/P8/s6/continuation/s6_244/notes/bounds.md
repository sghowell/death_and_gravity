# Whole all-momentum state and subtraction bounds

## 1. Mixed jets and the complete finite reference

For all admitted real histories and all k, set nu²=(99/100)(n+k²/16). Then nu<=omega<=6nu and nu>10^98. Use K0=10^90 only as a denominator majorant; there are NO omitted low modes.

Let delta=1/100. The noncommuting matrix-exponential derivative bound is exp(delta) times the product of its direction norms. Time/parameter Faà di Bruno therefore gives Bell majorants

B0=1, B_(j+1)=delta sum_(k=0..j) binom(j,k) B_k,
B_(j,a)=sum_(k=0..j) binom(j,k) a^(j-k) B_k, a=0,1,2.

These bound all mixed derivatives of exp(-Q-epsilon G), including noncommuting jets, with the common factor exp(delta). Products of a differentiated source factors supply a^(j-k). The case a0 uses0^0=1.

For a^-2=(1+t²)^-4, raw normalized derivatives are p_j(t)/(1+t²)^j, where p0=1 and p_(j+1)=(1+t²)p_j'-2(j+4)t p_j. Their full coefficient sums bound every real time in[-1,1]. The normalized frequency-square mixed jets have bounds

q_(j,a)=2 sum_k binom(j,k)c_(j-k)B_(k,a),

with q00=1 and the sharper q10=5. The factor2 exceeds exp(2delta), so it covers division by the full positive frequency square, including k0.

Differentiate omega² and omega*(1/omega)=1 with full mixed binomial coefficients. Their positive recurrences define U_(j,a) and V_(j,a), bounding |partial_t^j partial_epsilon^a omega|/omega and omega|partial_t^j partial_epsilon^a(1/omega)|. There are39 entries each through j12,a2. No derivative of a fixed denominator is omitted: the recurrence is the complete product rule.

The squeeze majorants S_(j,a) use lambda=omega'/omega and theta=3H+tr Q'/2. The background bound is |H|<=2, |H^(j)|<=4j! for j>=1. Source traces are bounded by three operator norms. Thus s0<5, |omega1|<=omega and |omega2|<=2omega.

The full coefficient-jet recurrence in estimates.py gives b_(n,j,a) with

|partial_t^j partial_epsilon^a r_n|<=b_(n,j,a) omega^-n, n+j<=11, a<=2.

It contains the derivative term, every quadratic partition and every mixed product weight. All195 coefficients used for n1..10 are recorded. Summing the entire finite reference and residual yields

|rhat_a|<=A_a nu^-1, |F_a|<=C_a nu^-10,

where A<(2.129,6.007,34.033), C<(3.432*10^16,1.273*10^18,5.050*10^19). Exact rational values, not these rounded displays, are used in every gate.

## 2. Entire actual error and its derivatives

Start with the full initial error <=10^90 nu^-11. Since the exact pure graph lies in the unit disk and rhat<1/100, the real part of the error generator is below6. On the unit time interval exp(6)<1000. Hence

|e0|<=1000[10^90/K0+C0]nu^-10 <10^40 nu^-10.

This bootstraps |r|<1/10. The exact tangent A_r=2iomega-2s r then has growth below exp(1)<3. The oscillatory 2iomega is not bounded by an exponentially growing modulus estimate.

The complete source-error equations are

e1'=A_r e1+[2iomega1-s1(r+rhat)-2s rhat1]e0-F1,
e2'=A_r e2+[2iomega2-s2(r+rhat)-4s1 rhat1-2s rhat2]e0
+[4iomega1-4s1 r-2s(r1+rhat1)]e1-F2.

The final bracket retains -2s e1². Initial derivatives e1,e2 are zero by the common germ. Define

M1=12+S01/K0+10A1/K0²,
M20=24+S02/K0+(10A2+4S01 A1)/K0²,
M21=24+S01/K0+20A1/K0².

Complete positive inequalities give

E1=3[E0 M1+C1/K0]<10^42,
E2=3[E0 M20/K0+10^42 M21+10*(10^42)²/K0^10+C2/K0²]<10^44,

with E0=10^40. Thus |e_a| is bounded respectively by10^40 nu^-10,10^42 nu^-9,10^44 nu^-8. The exact unrounded bounds are retained in the report.

## 3. Whole covariance and physical current

On |r|<=1/10 the complete covariance and its first three Frechet derivatives have bounds2,16,128,512 in operator norm with direction magnitude one. To verify these constants directly, write Sigma=N/(1-r rsharp). The denominator is at least99/100; reciprocal derivative bounds2,1,10,16 suffice. Full numerator derivative bounds2,4,2,0 give first/second/third covariance bounds10,32,158, respectively below16,128,512. The sharper row-sum bound ((1+.1)²/2+.1)/.99<2 gives the zeroth bound. This argument includes both off-diagonal entries.

Apply the whole first and second chain rules to Sigma(r)-Sigma(rhat). The report retains the resulting complete constants D0,D1,D2, including terms e0 rhat2, e1 rhat1, e1² and e0 rhat1². Their powers are nu^-10,nu^-9,nu^-8.

For the full spatial vertex, |delta_D omega²/omega²|<=2||D||F. Its first/second source derivatives are bounded by6 and30 times ||D||F from the complete quotient rule. Consequently the full N_D and its first two derivatives have bounds4,6,30. These include its trace/volume terms.

In J_D=-omega tr(N_D Sigma)/2, the half-trace cancels the harmless two-dimensional operator-norm factor. The complete three-factor product rule, |omega|<=6nu and all U0a give

|partial_epsilon^a[J_actual-J_reference]| <10^50 nu^(a-9)||D||F, a0,1,2.

All coefficients are recorded before rounding. The physical metric contact is present because BOTH the vertex and covariance have been differentiated.

## 4. Entire marker tail, not only its first omitted term

Introduce a marker zeta in the FINITE graph reference and use |zeta|=omega/Cstar with Cstar=10^6 independent of n. Take physical source derivatives at FIXED zeta before selecting this pointwise circle; no derivative of a moving contour is used.

The coefficient bounds in section1 scale with the pointwise omega, so the complete reference graph contour bounds are sums b_(j,0,a)/Cstar^j. They are below(2.129*10^-6,6.007*10^-6,3.404*10^-5). The full rational covariance has no pole on this disk. Complete current contour derivative bounds are below(10,100,1000)omega||D||F.

Only this diagonal current is even in zeta. Cauchy and the geometric tail therefore bound the ENTIRE reference current minus its orders0,2,4 by

2B_a Cstar^6 nu^-5 ||D||F, B=(10,100,1000).

Combine this with the exact-state/reference error, converting the powers nu^-9,nu^-8,nu^-7 to the weaker nu^-5 using nu>1. Every complete combined coefficient T_a satisfies2T_a<10^51.

## 5. Full radial integral and unchanged cutoff

With the full measure d³k/(2pi)³,

integral nu^-5 =64(100/99)^(5/2)/(6pi² n)<2/n.

Thus ALL current derivatives a0..2 of the actual-state-minus-adiabatic integral are below10^51/n times ||D||F and the corresponding source norms. After normalization this is below10^-940. It is not yet the full current: the complete finite action and fixed profile are added in local.md.

The auxiliary condition is exactly Omega_star<=K, K>=2sqrt(n), independent of history. At zero transfer the two memory legs have identical momentum condition; the single contact leg has that same condition. On its complement nu>=sqrt(99/100)K. The full radial Jacobian is at most4nu²dnu, since32(100/99)²/9<4. The remaining integral is therefore bounded by

2T_a/[(99/100)K²] <10^51/K².

The complete normalized cutoff error is below10^51/(kappa K²)<10^-748/K². Finite local/profile terms are full and unchanged on both sides. No finite-K pressure or current mean is set equal to its unprojected value.

Finally, smooth prepared G has zero jets at the initial endpoint. Repeated integration and Cauchy--Schwarz on the unit interval give ||G||C12<=||G||H13. Integrating ||D(t)||F uses ||D||L1<=||D||L2. These give the claimed zero-transfer time graph and the SAME graph for its cutoff tail, without pretending to control a nonzero-transfer kernel.
