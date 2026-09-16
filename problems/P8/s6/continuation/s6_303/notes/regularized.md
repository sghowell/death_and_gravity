# Exact integer-pole cancellation and a uniform massive-block bound

Work at mu=nu=1,25/4<=s<=16,-1<=t<=0 and0<ell<3.
Let cm0,cm1,bm0 be the complete small-channel S284 massive
triangle coefficient, its first dimensional derivative and the
massive bubble coefficient. Define the entire real small-channel
massive/rational block

S_m=cm0*(ell*J_t-Q_t)/2-cm1*J_t/2+bm0*(ell-L_t)
    +E_compact+(6ell+14)*T0-6*T1-ell*R0+R1.

The old finite Gram terms are not hidden in this definition: all
three are separately regular and bounded in notes/rate.md.
Every compact evanescent, LSZ and physical-pole term is included above.

## Cancel first, bound second

Writef(t)=t^2*S_m=A_J*J_t+A_Q*Q_t+A_L*L_t+A_0, where

A_J=t^2(ell*cm0-cm1)/2,
A_Q=-t^2*cm0/2,
A_L=-t^2*bm0,
A_0=t^2[ell*bm0+E_compact+(6ell+14)T0-6T1-ellR0+R1].

All four rational functions are regular atzero. The exact master jets
areJ0=1,J'0=1/6,Q0=L0=0,Q'0=L'0=-1/6.
Substitution givesf(0)=f'(0)=0 identically in s andell.
This cancels both apparent integer poles before any absolute estimate.

For the parameter denominatorA=1-t*x(1-x), putb=x(1-x).
On the entire stated real interval,1<=A<=5/4,0<=b<=1/4.
The differentiated master kernels are

J':b/A^2, J'':2b^2/A^3;
Q':-b(1-ln A)/A^2, Q'':b^2(2ln A-3)/A^3;
L':-b/A, L'':-b^2/A^2.

Together with the undifferentiated kernels these have magnitude at
most1. ForQ use0<=ln A<1/4. Integration preserves the bound.
The signs and zero jets are checked exactly.

## Whole rational derivative norms

For each A_j and its first two derivatives, factor the denominator
after exact cancellation. Only factors s,s+t-4,t-4,s-4 are allowed.
Their magnitudes on this domain are at least1; actual gaps are larger.
For the numerator polynomial use the coefficient l1 norm weighted by
|s|<=16,|t|<=1,|ell|<=3, divided by the exact numerical denominator.

The resulting bounds for(A,A',A'') are

A_J: (7972,14810,16784),
A_Q: (2422,4680,5548),
A_L: (60323/4,374714/3,552949),
A_0: (1108223764786/225,4558231312812/5,20687039523642988/225).

This is a whole-function algebraic norm, not a fit or grid estimate.
Unexpected denominator factors are rejected by the norm implementation.

The product rule boundsf'' by
sum_(j=J,Q,L)(|A_j''|+2|A_j'|+|A_j|)+|A_0''|.
Sincef(0)=f'(0)=0, Taylor's integral formula
f(t)=t^2 int_0^1(1-r)f''(rt)dr gives

|S_m|<=82748158895162527/1800 <4.6*10^13.

The endpoint is a continuous extension. Neither an individual pole nor
an individual Laurent contribution is bounded before its cancellation.
