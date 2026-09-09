# Actual subclock observers

For a future curve gamma(lambda) write u_dot>0 and
v=a*x_dot/u_dot, its physical spatial velocity. A mode covector
p=(-s*k/a,k*n), |n|=1, s=c or 1, contracts to

    p.gamma_dot=(u_dot*k/a)*(v.n-s).

The exact three-dimensional Gram identity proves
v.n<=|v|. Thus |v|<c makes BOTH contractions strictly negative.
On compact sampling support smooth strict inequalities give
a uniform separation. This uses the actual clock cone, not
only the physical-metric timelike condition |v|<1.

The global bound c^2>=1199/1215 is greater than (993/1000)^2.
Hence |v|<=99/100 gives contraction margin at least 3/1000,
apart from the displayed positive u_dot*k/a factor. The
class includes the comoving observer at all finite times.

For comparison, gamma_v(u)=(u,v*F(u),0,0) has physical metric
norm -1+v^2 and proper time sqrt(1-v^2)*u. Choose v=999/1000.
The positive-even-polynomial bound at |u|<=1/100 gives
c(u)^2<v^2 throughout that interval. There

    n=(c/v,sqrt(1-c^2/v^2),0)

is a real unit direction whose clock covector annihilates
the curve tangent. Its physical trajectory is nevertheless
timelike. This is a whole-interval example, not a point scan.

For subclock curves the signed wavefront set from S6.99
misses the conormal bundle. Pullback therefore exists, keeps
the time-frequency signs and preserves positive type.
The positive-type pullback result is
[Fewster, Theorem 2.2](https://arxiv.org/pdf/gr-qc/9910060);
its hypothesis is verified by the contraction above.

The interluminal example only fails this sufficient criterion.
It does not prove that every proposed restriction is impossible
or that an inequality is violated. That distinction is also
explicit in [Fewster--Pfeifer--Siemssen, Section 3.2](https://arxiv.org/pdf/1709.01760).
Their electromagnetic energy functional is not substituted for
our actual scalar observable.
