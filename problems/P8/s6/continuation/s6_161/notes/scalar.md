# Uniform boundary-value bounds for all full scalar loop shapes

Let v and z denote the light-pair channel and the complementary
heavy-route invariant. In the physical rectangle they lie
in [-2,6]. Put delta(x)=1-v x(1-x), so -1/2<=delta<=3/2.
For a bubble the MS finite kernel is
ell-integral_0^1 log(delta-i0) dx, ell=log(mF^2).

For 0<=v<=6, |delta|<=1. The integral of |log|delta|| is the
real part of the subtracted equal-mass bubble. For 0<v<4 it is
2-2r arctan(1/r), r=sqrt(4/v-1); for 4<v<=6 it is
2-beta log((1+beta)/(1-beta)), beta=sqrt(1-4/v).
It lies between zero and two, with continuous value two
at threshold. For -2<=v<=0 it is at most log(3/2)<1/2.
The logarithm's imaginary part is bounded by pi<4.
Thus the normalized bubble modulus is at most ell+6.

Group the two light Feynman parameters as (1-H)x,
(1-H)(1-x). For a triangle the remaining heavy parameter
is H. For a box split it as Hy,H(1-y). The full polynomials
and measures, with their Q factor suppressed, are

    triangle: Delta=delta(1-H)^2+MH,
              (1-H) dH dx / (Delta-i0);
    box:      Delta=delta(1-H)^2+MH-z H^2 y(1-y),
              H(1-H) dH dx dy / (Delta-i0)^2.

They follow from all four light-heavy external squared
masses being one. For repeated heavy routes z=0. No
large-M expansion is made in these formulas.

Write Delta=delta+(M-2delta)H+A H^2. In the box A is
in [-2,2]. For M>=32, on the extended interval [-1,1],
Delta'>=M-7>=3M/4. Delta(-1)<0 and Delta(1)>0, so
the unique zero in this interval can be used even when
it lies just outside the integration range. Change variable
q=Delta and put U=Delta(1). For triangles U=M; for boxes
M-3/2<=U<=M+1/2.

For triangles w(q)=(1-H)/Delta' obeys
|w(0)|<3/M and |w'|<3/M^2 on that interval. Subtract w(0)
before the boundary value:

    integral_delta^M w(q)/(q-i0) dq
      =w(0)[log(M/|delta|)+i pi 1_(delta<0)]
        +integral_delta^M (w(q)-w(0))/q dq.

The second term is below 4/M in modulus. Averaging x and
using the logarithm bound gives

    |Triangle| <(3 log M+22)/M.

For boxes put f=H(1-H), d=Delta', w=f/d. The endpoint
weights vanish. Integrating by parts at nonzero regulator
and then taking the boundary value gives integral w'(q)/q.
The exact derivatives are

    w'=f'/d^2-f d'/d^3,
    w''=f''/d^3-3 f' d'/d^4+3 f(d')^2/d^5.

On the extended interval, |f|<=2, |f'|<=3, |f''|=2,
|d'|<=4. At the worst allowed M=32 the resulting bounds
are |w'|<=160/(27 M^2)<6/M^2 and
|w''|<=704/(81 M^3)<9/M^3; they improve with M.
Subtracting w'(0), its regular remainder is below 10/M^2.
Use log U<log M+1 and the same averaged log/phase bounds:

    |Box| <(6 log M+52)/M^2.

For isolated delta=0 values these formulas mean limits.
The remaining logarithmic x singularity is integrable,
including at the double threshold root. No bound is
claimed on the divergent integral of |Delta|^-2 through
an unsplit real pole. Constant i0 can be retained during
the integration by parts: zero endpoint weights remove
the boundary terms and the subtracted numerator gives a
uniform integrable limit.

Finally Cmax=L+g/(M-6) and the complete scalar amplitude
obeys

    Bscalar=3/(2Q)[Cmax^2(ell+6)
      +4 Cmax g(3 log M+22)/M
      +4 g^2(6 log M+52)/M^2].

This remains a bound on the complex physical amplitude,
not just on its real part. The exact code uses the
rational cap M<2^1024, hence log M<1024, and the inherited
independent enclosure of ell.
