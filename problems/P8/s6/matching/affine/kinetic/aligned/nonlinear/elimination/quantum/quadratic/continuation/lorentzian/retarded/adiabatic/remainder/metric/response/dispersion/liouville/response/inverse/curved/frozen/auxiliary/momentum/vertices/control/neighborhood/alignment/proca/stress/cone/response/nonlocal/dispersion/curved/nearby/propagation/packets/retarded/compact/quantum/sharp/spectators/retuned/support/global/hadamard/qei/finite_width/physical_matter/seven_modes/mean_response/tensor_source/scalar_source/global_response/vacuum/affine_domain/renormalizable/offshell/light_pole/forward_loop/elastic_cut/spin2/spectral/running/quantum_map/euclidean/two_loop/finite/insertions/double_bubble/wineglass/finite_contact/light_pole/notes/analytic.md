# First-sheet continuation and complete radial estimates

Set Qnorm=16 pi^2, and work initially with real Euclidean external
momentum p. The scalar bubble subtracted at the inherited reference is

    I_R(Q) = -1/Qnorm integral_0^1 da log(1+a(1-a) Q^2).

For clarity below write w=a(1-a), so 0<=w<=1/4.
After combining one outer light and one heavy denominator,
shift q=k+xp for real p and use
Delta=1+x(M-1)+x(1-x)p^2 and Q=q+(1-x)p.
For the squared heavy denominator the same combination has
weight 2x and power three instead of power two.

## Continuation without an assumed complex contour shift

The rational/logarithmic expressions derived with real p are
holomorphic in the complex Hermitian ball ||p||^2<3+epsilon,
for a sufficiently small epsilon>0 and M>32. Indeed
Re Delta>=1+x(M-4-epsilon), and, for real q with y=q^2,
Re Q^2>=y/2-||p||^2. The logarithm argument stays in a
strict positive half-plane. After the mixed reference difference
is taken, the integrable majorants below are uniform on compact
subsets. The other displayed integrals already have those
majorants. They therefore define holomorphic functions of the
complex momentum components in a connected neighborhood.

The real-p expressions equal the Euclidean Feynman integrals
on a real open ball. Uniqueness of holomorphic continuation
extends those identities, including scalar dependence on p^2,
to this complex ball. This establishes the shifted representation;
it does not assume that an arbitrary complex momentum contour
can be moved without checking singularities.

Use the entire external embedding

    p(s)=(i(s+3)/(2 sqrt(3)), (s-3)/(2 sqrt(3)),0,0).

Literal checks give p(s)^2=-s and
||p(s)||^2=(|s|^2+9)/6<=3 on |s-1|<=2.
The strict gaps persist on a neighborhood of this closed disc.
The zero-external reference is the scalar value at p=0.
Although the chosen p(0) is a nonzero complex null vector,
the continued scalar amplitude there has the same invariant
value. A constant reference is in any event annihilated by OS.

## Anchored logarithm estimate

For r=(1-x)p, ||r||^2<=3 and |r^2|<=3. Thus

    Re[1+w(q+r)^2] >= 1+w(y/2-3)
                         >= (1+w y)/4,
    |(q+r)^2-y| <= 2 sqrt(3y)+3 <= 4 sqrt(y)+3.

The same lower bound holds along the straight segment to
1+w y. Integrating the logarithmic derivative along that
segment, and using w/(1+w y)<=1/(y+4), gives

    |log(1+w(q+r)^2)-log(1+w y)|
       <= D3(y)=(16 sqrt(y)+12)/(y+4) <= 7.

The final inequality follows from
7t^2-16t+16=7(t-8/7)^2+48/7>0.
In particular |I_R(q+r)| <= [log(1+y)+7]/Qnorm.
The decaying D3, not merely its constant upper bound, is
needed for the mixed reference difference.

## Independent elementary integrals

The exact radial measure is d^4q/(2pi)^4
=Qnorm^-1 y dy. Independent primitives in the source prove

    integral_0^infinity y/(y+1)^3 dy = 1/2,
    integral_0^infinity y^(3/2)/(y+1)^3 dy = 3pi/8.

For the latter set y=t^2 and use
3 atan(t)/4-5t/[4(1+t^2)]+t/[2(1+t^2)^2].
Both endpoints and derivatives are checked.
The inequality log(1+y)<=sqrt(y) follows by differentiating
t-log(1+t^2), whose derivative is (t-1)^2/(1+t^2)>=0.
Also log(y+2)<=log 2+log(1+y)<1+sqrt(y).
The elementary atan integral gives pi<4; its eight-term
alternating lower polynomial gives pi>3.

## Mixed sunset after its local reference

The six mixed refinements have coefficient -Lg and scalar
J_mix(s)=integral D_light^-1 H_Q^-1 I_R(Q).
It still has an overall local divergence. Subtract J_mix(0)
with the common regulator before removing it. After combining
the outer denominators, split the difference into

    [I_R(q+r)-I_R(q)]/(y+Delta)^2
    + I_R(q)[(y+Delta)^-2-(y+Delta0)^-2].

Here Delta0=1+x(M-1), Re Delta>=1 and
|Delta-Delta0|<=3x(1-x)<=3/4.
The first radial term is bounded by
integral (16 y^(3/2)+12y)/(y+1)^3 =6pi+6<30:
use y+4>=y+1 in D3.
The second is bounded by
(3/2) integral y log(1+y)/(y+1)^3 <=9pi/16<9/4.
Therefore |J_mix(s)-J_mix(0)|<33/Qnorm^2 over the
whole closed outer disc. This is a finite reference difference,
not a bound on the unregulated divergent J_mix alone.

## Same-heavy sunset

The three same-heavy refinements have coefficient g^2/2
and scalar integral D_light^-1 H_Q^-2 I_R(Q).
Its Feynman weight integrates as integral_0^1 2x dx=1.
Using Re Delta>=1 and the logarithm bound gives a radial
upper bound 3pi/8+7/2<5. Its unprojected magnitude is
therefore below 5g^2/(2 Qnorm^2).

## Nested remainder

The actual inherited inner OS kernel has
R(z)=-g/Qnorm integral_0^1 dx log(1+b(x)z)/z,
where 0<=b(x)<=1/M. On the unshifted real Euclidean
line z=y+1>0, monotonicity gives

    |R(z)| <= g/Qnorm log(1+z/M)/z.

There is no extra complex-strip factor two on this real axis.
For the outer heavy denominators, the entire embedding gives
Re[(k plus or minus p)^2+M] >= y/2+M-3 > (y+M)/2.
Both outer heavy channels and their symmetry factor give the
total coefficient g. Their complete unprojected bound is

    2g^2/Qnorm^2 integral_0^infinity
       y log(1+(y+1)/M)/[(y+1)^2(y+M)] dy
      < 4g^2/Qnorm^2.

For the final inequality drop M>=1, use log(y+2)<1+sqrt(y)
and the two anchored radial integrals: 1/2+3pi/8<2.
This deliberately conservative step retains full convergence.
The separate alpha Pi1_R part is already on-shell projected.
