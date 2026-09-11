# Full-profile derivatives and total variation

For the positive real root,

    s(x)=x(1+x^8)^(-1/8),
    s'=(1+x^8)^(-9/8),
    s''=-9x^7(1+x^8)^(-17/8),
    s'''=-9x^6(7-10x^8)(1+x^8)^(-25/8).

The full profile is used, not a small-time Taylor
approximation. It is odd, with limits -1 and +1,
and s' is nonnegative and at most one. Thus

    integral s' dx=2,
    integral (s')^3 dx<=2.

On each half line s' is monotone. Applying the
fundamental theorem to (s')^2/2 gives exactly

    integral |s's''| dx=1.

The only nonzero extrema of s'' have |x|^8=7/10.
There is one hump on each side, and s'' vanishes
at zero and both infinities. Hence

    integral |s'''| dx=4 max |s''|<=36,

using |s''|<=9, which follows by splitting
|x|<=1 and |x|>=1. These are analytic
total-variation statements, not grid bounds.

For x>=1 write s(x)=(1+x^-8)^(-1/8).
Convexity gives

    0<=1-s(x)<=1/(8x^8).

Oddness gives the other tail. In physical time,
the deviation from each asymptotic mass is
integrable. The derivative tails also make
both integration-by-parts boundary terms vanish.
