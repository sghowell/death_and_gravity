# Exact finite bounds and complete selected sum

b2=(n^3-6n^2logn+9n^2-6nlogn-9n-1)/(6(n-1)^5),
 b3=(n^4-12n^3logn+28n^3-36n^2logn-12nlogn-28n-1)/
     (12(n-1)^7),
 d=(n-1)b3-b2,
 chi_extra=g^2 E(n)/(16pi^2 kappa), E(n)=-32[d+2/n^4].

All primitives, endpoint formulas and equal-mass limits have exact
checks. The combined heavy-mass limit is n^2 E(n)->8/3.

For n>1, set r=z(1-z)/M. Pointwise0<(n-1)r<1 gives
b2/3<-d<b2. The bounds M<=n and r<(1-z)/(n-1) give

 1/(60n^2)<=b2<1/[6(n-1)^2].

Thus n^2>360 implies E>0 and E<32/[6(n-1)^2].
At the original parameters, exact rational checks imply

 0<chi_extra<g^2/(25kappa n^2),
 chi_extra/A0<n/(100kappa)<10^-604.

The S346 finite bound supplies
|chi_core|>g^4/(80pi^2 n^4). The exact original inequality

 160n^4/[6g^2 kappa(n-1)^2]<1

proves chi_extra<|chi_core|. Consequently

 chi_known_selected=chi_core+chi_extra<0,
 |chi_known_selected|/A0<|chi_core|/A0<10^-207.

This is a LOCAL analytic-origin coefficient in a fully fixed
comparison, not an above-threshold amplitude approximation.
