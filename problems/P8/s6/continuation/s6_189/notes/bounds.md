# Full joint-spacetime fourth-Sobolev estimate

All norms use dt dx and the tensor Frobenius inner
product. Define the standard joint norm

    ||h||H4^2=sum_(j+|alpha|<=4)
                 ||partial_t^j partial_x^alpha h||L2^2,

with each spatial multiindex counted once. No volume
weight is hidden in this norm; the a^3 factors have
already been included in the differential operator.

On I=[-1/2,1/2], the exact coefficients obey

    1<=a<=25/16, |H|<=8/5, |H'|<=4, |H''|<=12,
    0<R0<72, |R0'|<=120, 0<A<500000, |A'|<4.

For example H=4t/(1+t^2),
H'=4(1-t^2)/(1+t^2)^2 and
H''=-8t(3-t^2)/(1+t^2)^3.
The H bound follows on0<=t<=1/2 from
8/5-H=4(1-2t)(2-t)/(5(1+t^2)), and oddness handles
negative t. The remaining bounds follow from these
rational expressions and

    R0=24(1+7t^2)/(1+t^2)^2,
    R0'=48t(5-7t^2)/(1+t^2)^3,
    A'= -R0'/36.

The actual m=1000 gives the positive A lower bound;
R0<72 implies A>5m^2/12-2>0.

The three Laplacian terms have coefficient Euclidean
norm sqrt(3), while
Delta^2=sum_i partial_i^4+2sum_(i<j)partial_i^2 partial_j^2
has coefficient norm sqrt(15). Cauchy-Schwarz therefore
gives

    ||Delta partial_t^j h|| <=sqrt(3)||h||H4
                              <2||h||H4, j<=2,
    ||Delta^2 h|| <=sqrt(15)||h||H4 <4||h||H4.

This keeps all spatial derivatives and tensor components.
Applying the complete coefficient list of operator.md
gives

    ||Ddag D h|| <=(18817/125)||h||H4 <160||h||H4,
    ||L_A h|| <3900004||h||H4.

The first sum includes respectively
1,6|H|,4|H'|+11H^2,|H''|+7|HH'|+6|H|^3,
4,4|H| and4 for the seven displayed derivative groups.
The second is [1+3*(8/5)+2]*500000+4.

Using pi^2>9 and the fixed kappa=1e800,

    ||Qloc h||L2
      <[3900004+160/60]/(72 kappa)||h||H4
      <60000/kappa ||h||H4
      =6e-796 ||h||H4.

Compact smooth TT fields satisfy every integration
and adjoint hypothesis. By continuity the estimate
also holds in the corresponding closed Sobolev
subspace. It is not a pointwise stress bound or a
norm on arbitrary scalar/mixed metric modes.
