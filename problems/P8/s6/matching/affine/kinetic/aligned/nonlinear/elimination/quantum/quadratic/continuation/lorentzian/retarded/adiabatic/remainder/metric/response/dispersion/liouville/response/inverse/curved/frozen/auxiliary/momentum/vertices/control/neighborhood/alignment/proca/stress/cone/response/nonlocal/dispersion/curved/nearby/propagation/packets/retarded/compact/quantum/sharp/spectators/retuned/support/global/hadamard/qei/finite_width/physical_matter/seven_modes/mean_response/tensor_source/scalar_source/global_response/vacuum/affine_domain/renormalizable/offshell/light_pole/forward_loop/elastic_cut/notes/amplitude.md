# Continuous physical amplitude bounds and exact angular average

Let k=(s-4)/2, B=M+k and C=-lambda4+g/(M-s), with g=G^2.
The actual on-shell tree amplitude, independently rebuilt from S6.110, is

A0(s,z)=C+2gB/(B^2-k^2 z^2).

For 4<=s<=6, B>k>=0 and M>s. This is a real nonsingular amplitude.
It is nondecreasing in z^2 because its derivative is
2gB k^2/(B^2-k^2 z^2)^2. At fixed z its s derivative is

g/(M-s)^2
 -g[(1-z)/(M+k(1-z))^2+(1+z)/(M+k(1+z))^2]/2.

The bracket is at most 2/M^2. Therefore the derivative is strictly
positive, since (M-s)^(-2)-M^(-2)>0. These are continuous inequalities
over the whole angular/energy rectangle, not a sample grid.

It follows that the minimum on [4,6] is A0(4,0), the minimum on [5,6]
is A0(5,0), and the maximum is A0(6,1). Exact actual rational parameters
give

A0(4,0)>23lambda, A0(5,0)>42lambda, A0(6,1)<73lambda.

The entire physical amplitude is thus positive, and its square can be
bounded without subtracting nearly equal enormous terms.

For k>0 the exact even angular average Q=integral_0^1 A0(s,z)^2 dz is

Q=C^2+4Cg atanh(k/B)/k
  +2g^2/(B^2-k^2)+2g^2 atanh(k/B)/(Bk).

The code differentiates the full primitive, including both
J1=atanh(kz/B)/(Bk) and
J2=z/[2B^2(B^2-k^2z^2)]+atanh(kz/B)/(2B^3k),
and checks its zero anchor. Thus the formula includes the entire
nonlocal vertex, including interference terms. It does not assume
positivity term by term in this closed difference.

On this real branch 0<=kz/B<1. At threshold k=0, atanh(r)/r tends
to one, so Q continuously approaches (C+2g/M)^2. The threshold value
is handled explicitly rather than substituted into a zero-over-zero
expression. For actual parameters the closed expression is poorly
suited to floating evaluation; the rigorous enclosures use the positive
integrand and the exact amplitude extrema above.
