# Complete uniform all-momentum heavy-state stress proof

The estimates below apply to the entire exact selected state at the actual parameters. A finite-order local heat expansion alone is not the stress. Exact rational gates in estimates.py check every displayed constant.

## 1. Common complex domain and frequency

Use the complete real interval I=[-1,1], a(t)=(1+t²)², Amax=4 and Omega=sqrt(n+p²/16). All p>=0 are retained. On the union of complex radius1/64 disks about I, put v=1+t². Then Re v>=4095/4096>99/100, |v|<=8321/4096<33/16, and |Im v|/Re v<=130/4095<1/30. Therefore

Re(v^-4)>[(1-6/30²)/(1+1/30²)²]*(16/33)^4>1/32.

Thus Re omega²>=Omega²/2, |omega|<=8Omega, |omega|>=Omega/2, with the principal branch analytic. The argument of v^-4 has cosine greater than1/2. Since n is positive, the projection on that sector gives |p²v^-4/(n+p²v^-4)|<2, including p0 by continuity. Also |p²/a² / omega²|<=2, |a^-3|<2, |H|<5, |H'|<12. For lambda=omega'/omega=-Hz and U=3H'/2+9H²/4, |U|<100, |lambda|<10, |lambda'|<400. Consequently B0=-U-lambda'/2+lambda²/4 has modulus<400. These estimates depend on the analytic reference metric, not on analyticity of the full smooth QG1 scalar coefficient.

## 2. Entire finite WKB iterates and actual residual

Set w0=1 and

w(j+1)=sqrt(1+epsilon²/omega²*[B0-(log wj)''/2+lambda(log wj)'/2+((log wj)')²/4]).

Time domains lose r=1/1024 at each step. On the ball |w-1|<=1/100, |log w|<=1/50. Cauchy bounds make the bracket less than10^5. The auxiliary epsilon disk has radius R=Omega/10^8, so |epsilon²/omega²|<=4*10^-16: the map stays strictly inside that ball. For two inputs, log differences have norms at most2 times the input difference. One Cauchy loss bounds the bracket difference by(3/r²+10/r)*||w-v||<10^7||w-v||. The normalized map contraction is at most K|epsilon|²/Omega², K=10^8, and the first increment is at most the same quantity. The auxiliary disk contraction is at most10^-8. All square roots/logs use their branches near1.

At epsilon1, |wj-w(j-1)|<=q^j, q=K/Omega². Hence Wj=omega*wj are the full positive real iterates, not truncated high-order polynomials. For f6=exp(-i integral W6)/sqrt(2W6), the complete oscillator defect is r6=W7²-W6², so |r6|<=192K^7/Omega^12<10^60/Omega^12. Seven steps leave a complex width9/1024; the later ten integrations by parts use only10/4096 more width.

On the real interval W6>=.99omega>=.99Omega. Variation of constants in the normalized f6,conjugate(f6) basis preserves the exact Wronskian. On I its matrix norm is at most2*10^60/Omega^13, with interval length2. Since n>10^196, the full exponent is less than1/2 and exp(delta)-1<=2delta gives total Bogoliubov coefficient error below10^62/Omega^13. This includes the full exact residual, not only its first asymptotic coefficient.

## 3. Sampling and exact state selection

The sampling w=f² is the normalized compact bump from [the exact state definition](state.md), supported[-3/4,-1/2]. Z>1/216. With s=1-x² and complex radius s/8 in x, the change in1-x² has modulus at most17s/64<s/3, and Re[1/(1-z²)]>=3/(8s). Cauchy in physical time radius s/64 yields

||w^(j)||infinity<=216*86^j*(j!)².

All endpoint jets vanish; the maximum through j10 is below10^36.

Use an exact normalized mode S initialized to f6 at t=-1 only as a COMPARISON BASIS; the selected SLE is basis independent. Its physical c1 satisfies c1>=Omega/128 by the Wronskian and a³<=64. The phase-removed physical self-pair energy amplitude of f6 is

A=a^-3[(omega²-W6²)+d²+2iW6*d]/(4W6), d=3H/2+(log W6)'/2.

On the remaining complex domain, d has modulus<14. The actual q bound gives |omega²-W6²|<=384K; the huge actual mass makes this divided by W6 less than1. Thus |A|<100 independently of momentum.

Integrate integral w A exp(-2i integral W6) by parts ten times. Write the resulting amplitude as sum_r w^(r) A_(j,r). The recursion is A_(j+1,r)=D(A_(j,r)/(2iW6))+A_(j,r-1)/(2iW6). Since |1/(2W6)|<=2/Omega and each derivative loses1/4096, the sum of analytic coefficient bounds is at most100*8194^j/Omega^j. Therefore the full approximate c2 is below10^79/Omega^10. The exact-mode correction is below10^70/Omega^12; hence |c2|<10^80/Omega^10 and |c2|/c1<10^83/Omega^11<1/2.

For delta=sqrt(c1²-|c2|²), choose alpha=sqrt((c1+delta)/(2delta)) and beta=-c2/sqrt(2delta(c1+delta)). This is regular at c2=0, preserves CCR and minimizes the energy. Its magnitude obeys |beta|<=|c2|/c1: putting r=|c2|²/c1²<=1/4, the equivalent positive-square inequality follows from (1+2r)²(1-r)-1=r(3-4r²)>0 for r>0, with equality at0. The state is the actual SLE, not the comparison WKB mode or a finite-order adiabatic state. Fix the reference state once; transport its Cauchy covariance on later perturbations.

## 4. Five time derivatives and full momentum integral

For y=(sqrt(Omega)chi,chi'/sqrt(Omega)), Omega is fixed in time. The exact evolution matrix is [[0,Omega],[-(omega²-U)/Omega,0]]; the approximate matrix adds r6/Omega in its lower-left entry. The physical quadratic matrices are

B_rho=a^-3/(2Omega) [[omega²+9H²/4,-3H Omega/2],[-3H Omega/2,Omega²]],

B_P=a^-3/(2Omega) [[-n-p²/(3a²)+9H²/4,-3H Omega/2],[-3H Omega/2,Omega²]].

Use the vector1-norm and its induced matrix norm; all harmless two-dimensional factors are included below. On real times the normalized comparison vectors have norm below20. The full physical quadratic matrices have a sharp bound1000Omega on the complex reference domain. The exact and approximate evolution matrices and the physical rho/P quadratic matrices have time derivatives through order5 bounded by C*Omega, C=10^30. This follows from their complete analytic expressions and Cauchy radius1/4096; the potential includes the full r6 for the approximate mode. The residual derivatives through order5 are below10^85/Omega^12.

Use Y_j=20*j!*(2C)^j for normalized exact/approximate mode derivatives: ||y^(j)||<=Y_j Omega^j. For their difference use Z0=10^65 and

F_j=10^85 sum_r binomial(j,r)Y_(j-r),

Z_(j+1)=C sum_r binomial(j,r)Z_(j-r)+F_j.

Then ||(y_exact-y6)^(j)||<=Z_j Omega^(j-13). This follows from differentiating the full first-order forced system; it does not assume a small unweighted high-frequency propagator.

For a quadratic readout define Q_j(L,R)=C sum_(r+a+b=j) multinomial(j;r,a,b)L_a R_b. Exact-mode stress differences are bounded by2Q_j(Z,Y) Omega^(j-12). The SLE versus exact comparison is bounded by4*10^83 Q_j(Y,Y) Omega^(j-10). The largest coefficients through j5 have respectively251 and271 decimal digits, both below10^300. Since Omega>1, both contributions are bounded by their coefficients times Omega^-5.

The full W6 phase-cancelled readouts, with d=3H/2+(log W6)'/2, are

E(epsilon)=a^-3[W6²+omega²+epsilon²d²]/(4W6),

P(epsilon)=a^-3[W6²-omega²(1-2z/3)+epsilon²d²]/(4W6).

This auxiliary epsilon counts derivatives; at epsilon1 these are the full physical readouts. It is not inserted into a physical oscillatory phase. They are holomorphic on |epsilon|<=R. The full iterates give |W6|<=9Omega, |1/W6|<=3/Omega and |d|<14. For pressure retain |1-2z/3|<3, not the energy-only bound. Both complete readouts are bounded by(3/2)[81+192+196/10^16]Omega<10^4 Omega. Its coefficients through degree4 are exactly the required adiabatic0,2,4 terms; the phase cancels before this Cauchy argument. The even tail at epsilon1 is at most2*10^4*(10^8)^6/Omega^5. After the readout's one derivative of W6 there is still a separate time Cauchy margin1/2048; max_(j<=5) j!2048^j<10^20. Five time derivatives therefore add at most10^20, so the full tail is below10^74/Omega^5. No Cauchy bound is applied to an exponentially large complex phase.

The summed complete mode difference is therefore below10^301/Omega^5 for each rho/P derivative through5. The radial integral with the full measure is

(1/(2pi²)) integral_0^infinity p²(n+p²/16)^(-5/2)dp
=64/(6pi² n)<2/n.

Hence a deliberately looser complete state-dependent/subtraction bound is10^310/n. This is all momentum, with no ultraviolet integration cutoff.

## 5. Local heat stress and new profile budget

The independent arbitrary-spatial-dimension mode integral and full D-dependent covariant pole counteraction match exactly. At mu1, before64pi², the finite scalar local density is
(3/2-log n)n²+(log n-1)nR_old/3-2log n*a2,
a2=R_old²/72+(Riemann²-Ricci²)/180 up to the compact divergence.

The whole pole counteraction is varied in D+1 dimensions before the limit. Its finite evanescent metric terms must be included; component pole subtraction alone is not the same scheme.

On the complex domain, the simple-pole formula for H gives |H^(j)|<=4*j!*(64/63)^(j+1). In particular bounds5,6,12,30 suffice for j0..3. The actual mass satisfies394<log n<462, inherited as an exact real logarithmic interval from the fixed S239 parameters. The local mass-curvature and curvature-square rho/P polynomials are bounded before64pi² by40000n and10^6. Separate the constant vacuum term before taking time derivatives. Thus all local stress derivatives through5 are bounded by

[462n²+10^12(40000n+10^6)]/576<10^399.

Adding10^310/n gives a total below10^400, or B=10^-400 against kappa0. This is the heavy reference stress only, not a response operator norm.

For the new fixed heavy scalar retuning use [the exact fixed profile](profile.md). Rewriting it as its affine clock part plus V=1-T times the remainder avoids a large derivative bound for T. On the real clock tube |u|<=1,7/8<=X<=9/8, the affine part and its mixed four-jets are at most9B/8. A complex X radius1/64 gives remainder bounds6B*j!*64^j*(9/55)^1024<B/8. The separate mass-one light-vacuum V term is belowB/8 as well. Thus the complete profile four-jet bound is below2B. Its reference functions are only required to be smooth in u, with the proved first five derivatives.

The profile, the existing S239 finite extensions and the S238 classical changes can be combined under a NEW normalized four-jet budget10^-399. No older budget is silently reused. New quantum response and full interacting/bounce control remain open.
