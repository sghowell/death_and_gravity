# Complete current coefficient four-jets, not clock truncations

Put |u|<=10^-3, |N-1|<=10^-3 and X=N^-2. Use the exact S6.253
factorizations of the ENTIRE R/F differences, the unchanged S6.238
heavy additions, and the S6.250 full normalized heavy source. The
source-pinned S6.255 ancestry is replayed before this certificate.

For the old analytic switch, use a complex (u,N) bidisk of radius
10^-4 about each real center. Then |N-1|<=11/10000,

    |X-1| <= d(2+d)/(1-d)^2 < 1/400,  d=11/10000,
    |X| > 399/400,  |(1-X)/X| < 1/300,
    Re(X^2) > 199/200,  |X|^2 < 2,
    |(1+u^2)^-3| < 2.

Thus for T=X^1024/[X^1024+(1-X)^1024] and b=1024X^2 exp(-1024X^2),

    |1-T| < 2*300^-1024,
    |b| < 2048 exp(-1000) < 1.

The last inequality follows already from the positive quadratic term of
exp(1000). On this bidisk the exact rational tree scalar is bounded by
10^6 and the old canonical vacuum polynomial by 10^205. The factors
exp(-u^4) and 1-b are bounded by 2 each. The entire old R-clock
difference is below 8*300^-1024; the old analytic scalar-switch
difference is below 10^207*300^-1024. Cauchy's mixed-derivative factor
through total order four is at most 4!*10^16. Exact integer comparisons
give respective four-jet bounds 10^-2500 and 10^-2300. The switch
denominator is nonzero throughout the bidisk.
The code also bounds the two actual polynomial numerators coefficient
by coefficient, checks the full tree denominator 200(1+u^2)^12 and
the vacuum polynomial denominator 150, and uses the complex denominator
lower bound rather than a real-axis sample. All four derivatives of
X=N^-2 are independently enclosed on the whole real lapse strip.

Retain all three fixed vacuum constants. Their normalized total is

    5*1000^4/(128 pi^2 kappa)
    +(3/2-log(n))*n^2/(64 pi^2 kappa)
    +3/(128 pi^2 kappa).

Here kappa=10^800, 10^197<n<10^198 and log(n)<600, since the positive
cubic partial sum of exp(3) exceeds 10. The code replaces pi^2 by its
strict lower bound 9 and evaluates an absolute rational majorant below
10^-400. Even the gross overestimate 10^300 for this constant would
leave its switched four-jets below 10^-2200.

The two fixed total normalized profiles obey their unchanged C5 bound
2*10^-400. They are real smooth functions; NO complex analyticity is
assumed for them. Their entire contribution is
T[A+B(X-1)], A=-pressure, B=-(rho+pressure)/2.
The ordinary product rule, |X-1|<1/400 and the switch derivative bounds
give a four-jet majorant below 10^5 times the profile bound, hence below
10^-390. The elementary product-rule factor 16*2*(1+121)=3904 is already
below 10^5, since positive-order switch jets are below one here.
All profiles are retained as actual fixed coefficients, not varied as
live quantum means and not set to zero.

The S6.238 bound applies to all mixed (u,X) four-jets of the full heavy
additions on |u|<=1 and 7/8<=X<=9/8. Our whole real lapse strip lies
strictly inside that tube. The N derivatives of X=N^-2 through four
have bounds 3,7,25,121. The fourth Faà di Bruno coefficient sum is

    3^4 + 6*3^2*7 + 3*7^2 + 4*3*25 + 121 = 1027 < 10^4.

The lower-order and mixed sums are no larger. Hence the full heavy
four-jets in (u,N) are below 10^-2696. For the finite onepoint source,
retain -J1(1-T)/sqrt(kappa), J1=-g/(32 pi^2), g=1/8192. Its full four-jets
are below 10^20*300^-1024/10^400 < 10^-2900.

Combining the entire pieces conservatively proves

    |partial_u^i partial_N^j (R-Rclock)| < 10^-300,
    |partial_u^i partial_N^j (F-Ftree)| < 10^-300,
    |partial_u^i partial_N^j j_normalized| < 10^-2500,
    i+j<=4.

These are error intervals around the rational tree, not identities
erasing a finite correction. The full R is even in u, so R_u(0,N)=0
for every N in the strip; the literal full binding is differentiated
in the exact residual check.

Finally, apply outward rational interval arithmetic to every one of
the 15 mixed derivatives of each complete rational tree function
on the real strip. The exact maxima plus the nonzero error 10^-300
are strictly below 200 for R and 10^4 for F. The helper uses exact
integer/rational arithmetic and outward algebraic roots with adaptive
dyadic precision; proof calculations never use floating endpoints.
These bounds alone do not imply a background lifetime, quantum mean
equation or cutoff. The subsequent notes supply separate classical
continuation and full-matrix arguments.
