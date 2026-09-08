# Rank-two curl forms and exhaustive kinetic-matrix routing

All statements here concern the literal constant two-trace family in
the formulation. The full 8-component rolling trace mass Schur matrix
is separately derived from all connection equations. In its internal
two-dimensional field space,

    D=(3/8)[[1,-7],[-7,1]],
    D(1,-1)^T=3(1,-1)^T,
    D(1,1)^T=-(9/4)(1,1)^T.

It is invertible, with determinant -27/4, but is not positive definite.
This mass fact alone is not our propagating-kinetic verdict.

## Full-rank negative transverse forms

The two actual rolling background traces vanish. Their full first
variations involve only the scalar lapse, its time derivative and its
spatial gradient. There is no transverse trace source, metric vector
or scalar shift contribution. Hence the two transverse trace vectors
have velocity matrix Z/2 in the quadratic action. If symmetric Z is
invertible and not positive definite, it has a negative eigenvalue.
Both transverse trace directions propagate when Z is invertible;
the temporal equations and scalar/metric-vector constraints cannot
remove a transverse velocity direction in this decoupled block.
Thus such Z already fails the physical kinetic test.

This argument must not be reused indiscriminately for rank-one
zero-Schur forms: there a remaining connection multiplier can impose
the projected trace itself. Those charts have the separate complete
rank-one/zero-Schur derivation.

## Positive-definite Z still has a negative longitudinal direction

Use true exact-gradient changes to remove the lapse-time derivatives
from both trace mass sources. Write the longitudinal variables as
W_i^I=partial_i sigma^I and temporal fields as t^I. Setting lapse
velocity-induced shifts aside only while extracting the pure
sigma_dot submatrix, the exact extra velocity density is

    L = t^T D^-1 t/2 + q*(sigma_dot-t)^T Z*(sigma_dot-t)/2.

The temporal equation is
(D^-1+qZ)t=qZ*sigma_dot. Where regular, elimination gives the exact
configuration-velocity coefficient

    K_sigma = (D+Z^-1/q)^-1/2.

The code derives it by substituting the temporal solution into the
unreduced quadratic form, then independently checks this inverse
expression. It also checks the determinant identity relating the two
temporal/elimination charts, so one is not inverted while the other
silently has a null mode.

Let m=(1,1)^T and p=(1,-1)^T. Positive-definite Z has a positive-definite
inverse. For the explicit sufficient bound

    q > (2/9)*m^T Z^-1 m,

we have m^T(D+Z^-1/q)m<0, whereas
p^T(D+Z^-1/q)p=6+p^T Z^-1 p/q>0. A real symmetric 2-by-2 matrix
with one strictly positive and one strictly negative quadratic value
has one eigenvalue of each sign and is invertible. Its inverse, and
therefore K_sigma, have the same inertia. No high-q series remainder
or extrapolation from a grid is needed.

The full scalar shift equation remains the original
v_dot=Theta*n-ell*s/2: neither exact-gradient vector sector adds a
scalar shift source. At every finite u!=0 this fixes n. After the
joint temporal/lapse/shift elimination, the restriction of the full
four-scalar velocity form to v_dot=s_dot=0 is exactly K_sigma.
The presence of a negative value in this restriction proves a
negative value in the full physical configuration-velocity form,
regardless of the remaining clock/matter mixing. This uses a regular
punctured chart, not a momentum/coordinate exchange at the center.

There are no leftover scalar multipliers hiding this negative direction:
the full (n,shift,t1,t2) auxiliary Hessian has determinant
-4*q^2*Theta^2*det(D^-1+qZ), independent of its lapse diagonal and
lapse/temporal cross entries. The exact block identity and two
independent joint solves on the actual positive- and negative-time
backgrounds are checked. All four auxiliaries are regular in the
strict domain above.

For Z=[[2,1],[1,3]], the sufficient threshold is 2/15. Taking q=2
provides a nonunit off-diagonal exact test with an indefinite
longitudinal velocity matrix. This is a momentum-domain example,
not an eigenfrequency or an EFT cutoff.

## Exhaustion of symmetric 2-by-2 constant forms

Write Z=[[r,t],[t,s]]. The zero matrix is the auxiliary control.
If det(Z)!=0, the two full-rank cases above exhaust its inertia.
If det(Z)=0 but Z!=0, then:

- for r!=0, Z=r*(1,t/r)^T*(1,t/r);
- for r=0, det(Z)=0 implies t=0, and
  Z=s*(0,1)^T*(0,1) with s!=0.

These are exactly the rank-one family, including both signs of its
single nonzero kinetic eigenvalue and its gamma=0 Schur charts.
No nonzero symmetric constant matrix is left unassigned. This
finite rank routing does not classify other connection tensors,
field-dependent coefficients or additional operators.
