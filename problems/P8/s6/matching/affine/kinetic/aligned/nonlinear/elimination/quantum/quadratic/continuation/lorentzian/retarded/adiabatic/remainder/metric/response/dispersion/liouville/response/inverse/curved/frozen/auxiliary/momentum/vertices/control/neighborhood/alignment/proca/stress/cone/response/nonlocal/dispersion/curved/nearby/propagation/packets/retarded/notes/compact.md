# Smooth compact local probes with matter-spacelike separated supports

Choose interior source and detector times

    s0=-T/4, t0=T/4, T=10^-7,

and x0=(S_c(t0,s0),0,0), y0=0. The actual whole-interval
speed bounds from S6.89 give

    S_c(t0,s0)-S_m(t0,s0)>d=99/(8*10^14).

Let h=d/100. Take source and detector neighborhoods with time
half-width h and Euclidean spatial-ball radius h about
(s0,y0) and (t0,x0), respectively. They stay strictly inside
I. For every source/detector pair in them, the spatial
distance can decrease by at most 2h. The matter radius can
increase by at most 4h because omega_m<2 and each time
endpoint changes by at most h. Thus

    |x-y|-S_m(t,s)>d-6h>0, t-s>(t0-s0)-2h>0.

support.py verifies the exact margins, as well as two other
interior time pairs. In this flat-spatial homogeneous physical
metric the matter null radius is exactly S_m, so every pair
of points in the two neighborhoods is matter-spacelike.
This is stronger than comparing only their centers.

The clock-front distribution passes through the interior
center pair. Its positive surface delta cannot be cancelled
by the locally L2 remainder, and the matter front is absent
throughout these neighborhoods. Consequently the exact
retarded kernel restricted to their product is nonzero.

There therefore exist real smooth compactly supported J in
the source neighborhood and f in the detector neighborhood
with nonzero retarded pairing integral f G J. Here is an
elementary separation argument: if every such product test
had zero pairing, convolve the kernel with a product of
source and detector compact smooth approximate identities.
Every resulting smooth value on a smaller product
neighborhood would vanish. Taking the approximate-identity
limit would make the distribution itself vanish there,
contradicting the nonzero clock shell.

These are genuinely compact local test supports, unlike the
S6.90 packets with retained tails. The statement is existence
of some probes, not a specified pulse with a certified
finite-band preparation, magnitude, noise or detector
resolution. An arbitrarily small rescaling of J keeps a
nonzero linear response, but does not prove nonlinear
control of a physical source apparatus.

Physical-volume smearing factors are positive smooth
functions and may be absorbed into f without changing
support. The result is independent of that smearing
convention. It concerns a linear classical relational
observable on the fixed actual background, not a
nonperturbative local observable of quantum gravity.
