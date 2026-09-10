# Exact whole-tube mixed C4 error bounds

Use the same external h and h^2 weights, normalized mixed
derivatives and full retuned lower function as S6.108.
Only the switch majorant changes. On 9/10<=X<=11/10 set

    q=(1-X)/X, |q|<=1/9,
    |q^(j)|/j! <=(10/9)^(j+1), j>=1.

For j<=6 the normalized q^n derivative bound is the
coefficient of z^j in the finite Taylor expression

    sum_(l=0..j) binomial(n,l) (1/9)^(n-l)
       [sum_(i=1..6)(10/9)^(i+1) z^i]^l.

Call it D_j. An independent rational differentiation
recurrence verifies every coefficient formula through order
six; it does not sample the high-degree switch.

Because n is even, 1+q^n>=1 pointwise. The finite Taylor
inverse jet norm for its reciprocal is bounded by the
finite geometric sum in the nonconstant D jets. The native
sum_(j<=6) D_j<1/2 implies an inverse norm below 2.
No expansion of enormous powers of rational denominators,
nor an assumed infinite-series convergence, is needed.

For w=nX^2 exp(-nX^2), strip the exponential:
P_0=nX^2, P_(j+1)=P_j'-2nX P_j. Bound each polynomial
monomial at |X|<=11/10. The positive Taylor terms show
exp(81/100)>2, hence exp(-nX^2)<2^-n on the tube.
The resulting six-jet norm of w is below one. Therefore
each normalized S derivative through order six is bounded by

    4 sum_(i=0..j) D_i.

Use S4=sum_(j<=4) S_j and Sd4=sum_(j<=4)(j+1)S_(j+1).
The same exact difference identities for F2,A3,A4,A5 then
apply with E=S+(X-1)S' and F=(X-1)S. The inverse-h weighted
norm is 210, both tube R factors exceed 4/5, and the new R
jet norm is below 233. The unchanged finite inverse-R bound
and the positive scalar envelope n L+C from S6.108 give
all five weighted mixed C4 coefficient errors below 10^-400.

Increasing n by two contracts every D jet, including the
extra affine scalar n factor, by at most
(1/81)(1+2/(1024-5))^7<1/50. The bump majorants contract by
at most (1/4)(1+2/1024)^7<1/3. Thus the same ceilings are
available for every even n>=1024. These are exact rational
coefficient bounds on the full time tube, not grids or
quantum error estimates.
