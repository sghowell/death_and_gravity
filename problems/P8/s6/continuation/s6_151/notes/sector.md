# Equal-mass derivative and endpoint-resolved parameter sectors

Let e=epsilon, d=4-2e, Q=16pi^2 and all three squared masses
equal b. The Euclidean sunset integral has Schwinger form

    S(b)=exp(2 gamma e) mu^(4e) Gamma(-1+2e)
         b^(1-2e) A(e)/Q^2,
    A(e)=integral_simplex (xy+xz+yz)^(-2+e) dx dy.

Differentiating the common squared mass inserts a squared
propagator in each of the three symmetric lines. Thus the
zero-momentum wineglass is -S'(1)/3, giving

    J_bare(0)=exp(2 gamma e) mu^(4e) Gamma(2e) A(e)/(3Q^2).

The derivative is justified on the common absolute-convergence
domain 1/2<Re e<1, where both massive integrals and their
derivatives converge. The identity then continues meromorphically.

Divide the simplex into three sectors according to its largest
coordinate; divide each again according to the other two.
In one sector set the projective coordinates to (1,r,rt),
0<=r,t<=1. The six equivalent sectors give

    A(e)=6 integral_0^1 dt integral_0^1 dr r^(e-1) B(r,t,e),
    B=[1+r(1+t)]^(1-2e)/[1+t+rt]^(2-e).

Subtract b0(t,e)=(1+t)^(-2+e) at r=0. Then

    A(e)=6[P(e)/e+R(e)],
    P(e)=(1-2^(e-1))/(1-e),
    R(e)=integral dt dr r^(e-1)[B-b0].

The apparent P singularity at e=1 is removable. On compact
epsilon sets near zero, B-b0=O(r) uniformly in t, and its
epsilon derivatives add only bounded logs of denominators in
[1,3]. Factors log r from r^e remain integrable for Re e>-1.
R is therefore holomorphic there and can be differentiated
under the integral. This resolves the actual endpoint poles,
rather than assigning a number to a divergent integral.

At e=0 write v=1+t+rt, w=1+r(1+t), B00=w/v^2,
b00=(1+t)^-2 and D=(B00-b00)/r. Exact radial integration,
followed by u=t/(1+t), turns R(0) into the integral from
0 to 1/2 of

    -1+1/[2(1-u)]+3/[2(1+u)]-log(1+u).

A primitive is 1-log(1-u)/2+(1/2-u)log(1+u).
It gives R(0)=log(2)/2. Since P'(0)=(1-log2)/2,

    A(e)=3/e+3+6 j e+O(e^2),
    j=(1-log2-(log2)^2/2)/2+R'(0),

where

    R'(0)=integral dt dr {
      log(r) D+[B00(log v-2log w)-b00 log(1+t)]/r }.

This is a convergent expression for the finite constant.
Its diagnostic value is about -1.49442065276884266.
The bound, not that decimal, is used in the certificate.

Independently, integrate the inner bubble's Feynman parameter
and the outer radial variable y. Setting y=z/(1-z) gives
Euler weights z^(1-e)(1-z)^(2e-1). The result is

    A(e)=3/[e(1+e)] integral_0^1
      2F1(e,2-e;2+e;1-x(1-x)) dx.

This second representation checks the normalization and all
six sectors without using their subtraction algebra.
