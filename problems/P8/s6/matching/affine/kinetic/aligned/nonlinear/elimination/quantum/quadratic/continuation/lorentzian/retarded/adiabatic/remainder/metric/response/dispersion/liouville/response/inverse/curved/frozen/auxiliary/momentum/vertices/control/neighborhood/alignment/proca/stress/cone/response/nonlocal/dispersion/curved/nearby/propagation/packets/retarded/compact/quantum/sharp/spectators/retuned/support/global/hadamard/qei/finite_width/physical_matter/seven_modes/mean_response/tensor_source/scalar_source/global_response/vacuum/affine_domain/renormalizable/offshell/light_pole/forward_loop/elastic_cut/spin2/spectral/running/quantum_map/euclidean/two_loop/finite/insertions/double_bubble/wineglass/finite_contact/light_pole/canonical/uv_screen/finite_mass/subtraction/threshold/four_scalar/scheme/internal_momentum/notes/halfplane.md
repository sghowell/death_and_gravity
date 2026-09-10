# Complete subtraction and an unbounded complex-half-plane bound

Write Q=16 pi^2, N=6, Y=y^2, A=x(1-x), and m=mF. The inherited full finite
fermion kernel at its named scale is

    f(s) = (2 N Y/Q) [2m^2 + integral_0^1 g_x(s) dx],
    g_x(s) = -(4m^2-s) log(1-A s/m^2).

The complete dimensional pole is affine in s. Subtract f(1) and
(s-1)f'(1), including the finite parts. Every local affine term cancels
before an estimate. The huge finite mass reference is not used as a
negative free scalar mass. All internal free scalar lines retain mass one.

Independently differentiating the full kernel gives the useful identity

    g_x''(s) = -A/(m^2-A s) - m^2 A(1-4A)/(m^2-A s)^2.

For 0<=x<=1 both numerator weights are nonnegative, with integrals 1/6
and 1/30. Fix 1<=sigma<4m^2 and d_sigma=m^2-sigma/4>0. On
Re(s)<=sigma, Re(m^2-A s)>=d_sigma. The principal logarithm is analytic
there. The segment from one to s stays in this half-plane, including
arbitrarily large imaginary part and negative real part.

Taylor's exact integral remainder therefore gives

    |f_R(s)/(s-1)^2|
      <= (N Y/Q) [1/(6d_sigma)+m^2/(30d_sigma^2)] = B_sigma.

The quotient extends removably at s=1. This is not a truncated expansion
in internal momentum divided by m. The finite gap and integrable
parameter weights justify differentiation and integration on each compact
subset; the displayed majorant is uniform on the entire half-plane.

At the actual parameters, with sigma=1, both the unscaled B_sigma and its
division by the full canonical kinetic lower bound are below 10^-606.
The bound is sharper than the inherited small-disc coefficient. The latter
is not changed or invalidated.

The unscaled coefficient is the formal one-loop insertion. Dividing by
the fixed one-loop kinetic normalization refers to the selected quadratic
functional, not an all-orders resummation. A constant uniform insertion
bound alone does not establish convergence or the value of an outer loop.
