# Full physical-angle and forward bounds

Let s0=4/3 and epsilon_vertex=8(4-n)/(16pi² k). The complete new contribution after the new value condition is exactly

A_extra,OS4=epsilon_vertex*g²[sum_channels1/(n-channel)-3/(n-s0)]

-[K2(sigma2-16/3)+K3(sigma3-64/27)]/(16pi²).

No heavy-mass expansion is used in this expression. Combining it with the already controlled complete polynomial loop is justified by the exhaustive graph classification, not by identifying the two actions.

## Finite contact and contact-only potential

For0<=x<=1, n>2, 1<=F(x)=(1-x)²+nx<=n. The independent mixed representation gives0<=Kmix<=3 log n<1386. The full new symmetric value is bounded by

B_extra=[32g²*1386/k+64g²/(k(n-4/3))+|K0|+(16/3)|K2|+(64/27)|K3|]/144,

using pi²>9. This bounds the complete diagrams, not an endpoint expansion. The old negative symmetric value has positive magnitude between5985g^4/(640D²) and14793g^4/(192D²). The exact new bound is smaller than that lower bound; hence

0<deltaC_new<14793g^4/(192D²)+B_extra<10^-408.

Its extra contribution divided by24q is below10^-200. Subtracting this finite contact from S238's completed-square quartic leaves beta-deltaC_new/24>q/2. Thus the **contact-only** comparison still has V_rel>=H²/6+qPhi^4/2. This does not assert a minimum of the full bare potential after the finite linear and quadratic terms, or a quantum effective-potential theorem.

## Complete physical window

For4<=S<=10^196 and all physical angles, the centered exchange is nonnegative. With a=n-s0 and x_i=channel_i-s0, use the exact identity

1/(a-x)-1/a-x/a²=x²/[a²(a-x)]

and sum_i x_i=0. The centered exchange times g² is the full separate tree minus its strictly positive symmetric tree value4g²/[D²(3D+2)]. It is therefore at most that tree, itself at most61/60 of the original tree. We use the weaker factor2.

Also |sigma2-16/3|<=4S², |sigma3-64/27|<=S³, and A_original>=lambda S²/4. All denominators remain positive since S<n. Therefore

|A_extra,OS4|/A_original
<=2*[8(n-4)/(144k)]+[16|K2|+4|K3|*10^196]/(144lambda)
<10^-600.

Adding the old complete first-loop bound10^7g²/(9n) still gives a strict total below10^-199. With the old separate-tree error, the complete tree-plus-first-loop error is below1/60+10^-199<1/59. This is explicitly a truncated perturbative matching result, not an omitted-loop estimate.

## Forward coefficients

Use v=s+t/2-2, u=2-v-t/2 and define b20=[v²]A, b21=[v²t]A, b40=[v^4]A, with derivative divisors2,2,24 respectively. The exact new shifts are

delta_b20_extra=epsilon_vertex*(4lambda)-2K2/(16pi²),

delta_b21_extra=epsilon_vertex*(-3gamma)+K3/(16pi²),

delta_b40_extra=epsilon_vertex*(gamma²/lambda).

In particular the new b40 shift is **not zero**. The full rational exchange is analytic near the subtraction region because its heavy poles are outside it; these derivatives do not follow merely from a real-axis error bound.

Absolute relative bounds are respectively epsilon_upper+|K2|/(288lambda), epsilon_upper+|K3|/(432gamma), epsilon_upper, with epsilon_upper=8(n-4)/(144k). Each is below10^-600. The old strict b20 and b21 lower bounds dominate these errors. Thus the complete new first-loop shifts satisfy

0<delta_b20/(4lambda)<10^-203,

0<delta_b21/(-3gamma)<10^-203,

|delta_b40|/(gamma²/lambda)<10^-192.

The sign of the total b40 shift is unassigned; the positive tree-plus-first-loop b40 margin remains. All actual comparisons use exact rationals, not floating-point underflow.
