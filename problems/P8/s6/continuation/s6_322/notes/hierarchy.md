# Complete temporal pair: hierarchy and compatible faces

Let u=a+b, W=u+c, Q_A=q1+q2 and q_c=c*n_c, with fixed real null
directions and unit complex TT leaf tensors. The complete temporal pair
is H2 and F2=ab*H2/u. The all-angle exact coefficient calculation gives

 ||F2||<=530*u/sqrt(kappa),
 ||partial_a F2||<=1566/sqrt(kappa),
 ||partial_b F2||<=586/sqrt(kappa),
 ||partial_ab F2||<=3104/(u*sqrt(kappa)).

For z=a/u, each component of F2 is u*P(z,r)/sqrt(kappa), where the
angular denominator is a positive constant times(1+r^2)^m and each
angular numerator degree is<=2m. Differentiation gives
 P+(1-z)P_z, P-zP_z and -z(1-z)P_zz/u.
Absolute coefficient sums on0<=z<=1,r>=0 prove the bounds. All24
component identities, their axis limits, and equality to the complete
temporal current were checked exactly in probe27536 (C0,4.331247292s).
Unit complex TT leaves follow by the same l1 coefficient argument as
S317. The compatible origin is zero by the O(u) bound.


The six independent coefficients are(T11,2T12,2T13,T22,2T23,T33).
Thus their absolute sum bounds the full matrix entrywise l1 norm,
including both off-diagonal entries, and hence the Frobenius norm.
This identity is replayed against the actual complete temporal pair,
not inferred only from contraction against one chosen source.

The normalized expressions are
P, P+(1-z)P_z, P-zP_z, -z(1-z)P_zz.
There are24 coefficients,96 derivative coefficient certificates and
48 exact single-axis limit identities. The O(u) norm bound supplies
joint-origin compatibility, including arbitrary energy hierarchies.
Exact collinear original poles remain outside the point evaluator.

