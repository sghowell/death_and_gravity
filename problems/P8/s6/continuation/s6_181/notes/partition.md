# Compact dyadic bands and their Fourier L1 bound

Choose chi(y)=1 for y<=1, chi(y)=0 for y>=2, and on 1<y<2
use chi(y)=1-10t^3+15t^4-6t^5 with t=y-1. Its first two
derivatives vanish at both transition endpoints, so chi is C2.

The exact derivative -30t^2(1-t)^2 shows monotonicity and
|chi'|<=15/8. Write u=1-2t. Then chi''=-15u(1-u^2).
For R=u^2 in [0,1],

    4/27-R(1-R)^2=(R-1/3)^2(4/3-R)>=0,

so |chi''|^2<=100/3<36 and |chi''|<6.

For every integer j put h=2^j and

    phi_j(x)=chi(x/h)-chi(2x/h).

Monotonicity gives 0<=phi_j<=1. Its support is [h/2,2h]
of length 3h/2. The two transition derivatives have disjoint
interiors: one is on [h/2,h], the other on [h,2h].
Therefore the derivative bounds are maxima, not an unjustified
cancellation:

    |phi_j'|<=15/(4h), |phi_j''|<=24/h^2.

Finite sums telescope to chi(x/2^N)-chi(2^(N+1)x);
for every fixed x>0 this is exactly1 once N is large.
Thus the full all-integer partition sums to1 and is locally finite.

Suppose max(|a|,h|a'|,h^2|a''|)<=A_j on the band.
Then g_j=a phi_j is compactly supported C2 with g_j and g_j'
zero at the two outer endpoints. Product bounds give

    G0=||g_j||L1<=3h A_j/2,
    G2=||g_j''||L1
       <=(3/(2h))[1+2*(15/4)+24]A_j
       =195 A_j/(4h).

Twice integrating by parts in its Fourier transform gives
|g_hat_j(t)|<=min(G0,G2/t^2). This is valid for t>0;
there is no claim of an absolutely convergent unweighted
full-frequency integral. The band sine kernel is

    K_j(t)=-2 Im[e^(2it) g_hat_j(t)].

Splitting the time integral at sqrt(G2/G0) (zero bands are
trivial) gives

    ||K_j||L1(0,infinity)<=4sqrt(G0 G2)
        <=sqrt(1170) A_j<35 A_j.

Both the sine factor two and the half-line time domain are
retained. No distributional delta function is inserted.
