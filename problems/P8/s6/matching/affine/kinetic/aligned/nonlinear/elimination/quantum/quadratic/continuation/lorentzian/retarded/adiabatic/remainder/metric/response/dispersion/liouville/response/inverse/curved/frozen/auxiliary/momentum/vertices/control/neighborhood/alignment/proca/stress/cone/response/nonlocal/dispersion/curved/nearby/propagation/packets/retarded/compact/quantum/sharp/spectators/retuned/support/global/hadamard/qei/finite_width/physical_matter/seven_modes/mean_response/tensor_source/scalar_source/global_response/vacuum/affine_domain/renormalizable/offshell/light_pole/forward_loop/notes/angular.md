# Uniform angular remainder over the entire loop integral

Write c=g/M^2, d=lambda4-3g/M. The actual parameters obey 0<d<5c.
For a channel z=s or u, C(z)=-lambda4+g/(M-z). After the light-pair
shift the vertex is V=C+g[(A+2q.r1)^(-1)+(A+2q.r2)^(-1)],
A=y+M-ell. The t-channel has ell=1, C=C(0), and external shifts
+/-p1 at one vertex and +/-p2 at the other.

On the declared disc, |A|>=(M+y)/2 and ||r_i||<=sqrt(3/2).
For z_i=2q.r_i,

|z_i/A| <=4 sqrt(3/2) sqrt(y)/(M+y)
 <=sqrt(6/M)<1/2,

uniformly for **all** y>=0. The last inequality follows by maximizing
y/(M+y)^2, whose maximum is 1/(4M). Use the exact identity

1/(A+z_i)=1/A-z_i/A^2+z_i^2/[A^2(A+z_i)].

This expands only the external shift; it never expands y/M.
Decompose V=V0+L+Q, with V0=C+2g/A radial and L odd in q.
Conservative uniform majorants are

|L| <=20 g sqrt(y)/(M+y)^2,
|Q| <=192 g y/(M+y)^3.

The linear constant bounds 16 sqrt(3/2)<20. For each Q term,
|z_i|^2<=6y, |A|^2>=(M+y)^2/4 and |A+z_i|>=(M+y)/4;
summing two terms gives 192. The ratio of these Q and L majorants
is at most 24/(5 sqrt(M))<1 because M>24.

The small radial vertex must retain the low-momentum cancellation.
Exactly,

V0=-d+g s/[M(M-s)]-2g(y-ell)/(MA).

The first two terms have absolute value below 11c. The last is bounded by
7c M(y+1)/(M+y). Since M(y+1)/(M+y)>=1, the common bound

|V0| <=20c M(y+1)/(M+y)

is valid. Merely bounding each original coupling separately would lose
this crucial cancellation.

Average over the unit three-sphere of real shifted q. The linear
V0 L terms vanish by q->-q. For two possibly different vertices with
the same bounds (in particular the t channel), the remainder of their
product relative to V0_left V0_right obeys

|average(V_left V_right)-V0_left V0_right|
 <=2 V0_bound Q_bound+(L_bound+Q_bound)^2
 <=7680 c^2 M^3 y(y+1)/(M+y)^4
    +1600 c^2 M^4 y/(M+y)^4.

No complex-conjugation assumption is used. Divide by the light denominator
absolute square, bounded below by (y+1/4)^2, and include radial measure y.
Positive gap polynomials prove

y^2(y+1)/(y+1/4)^2 <=4y,
y^2/(y+1/4)^2 <=1.

The exact convergent radial integrals are
integral y/(M+y)^4 dy=1/(6M^2) and
integral 1/(M+y)^4 dy=1/(3M^3), over zero to infinity.
Thus one channel before its loop prefactor is bounded by
(5120+1600/3)c^2 M. Summing all three channels and using pi>3 gives

sup_disc |angular remainder amplitude|
 < (530/9)c^2 M <60c^2 M.

This remainder is holomorphic by the complete parameter-domain argument
and its uniform integrable majorant. Cauchy's coefficient estimate on
the unit-radius circle gives the **same** upper bound on its b2.
The finite local subtractions are radial and do not alter this angular
remainder. The t-channel angular dependence is included, not dropped.
