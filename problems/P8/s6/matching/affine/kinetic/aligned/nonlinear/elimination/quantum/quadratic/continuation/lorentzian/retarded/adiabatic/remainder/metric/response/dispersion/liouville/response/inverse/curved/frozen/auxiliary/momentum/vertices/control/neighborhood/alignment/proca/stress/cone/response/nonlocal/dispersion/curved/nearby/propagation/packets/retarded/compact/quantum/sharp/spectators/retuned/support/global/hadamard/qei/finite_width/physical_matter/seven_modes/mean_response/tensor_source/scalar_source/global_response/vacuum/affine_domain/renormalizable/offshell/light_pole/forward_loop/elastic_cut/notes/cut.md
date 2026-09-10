# Strict integrated weight and the cut-subtracted coefficient

The independent arctangent identity pi/4=integral_0^1 dt/(1+t^2)
and positive gap t^2/(1+t^2) prove pi<4 strictly. The parent's
independent positive-integral argument proves pi>3.

On 5<=s<=6, beta>=2/5 and A0>42lambda.
On all 4<=s<=6, beta<=3/5 and A0<73lambda.
Using the actual squared angular integral and 3<pi<4,

rho(s)>441 lambda^2/80>5lambda^2 on [5,6],
0<=rho(s)<5329 lambda^2/160<34lambda^2 on [4,6].

At threshold rho=0; away from it rho>0. The strict lower bound uses
a positive-length interior subwindow, not a single sample point.

For the forward Taylor coefficient b2 (half the second derivative),
the low-energy cut functional is
I_[4,6]=(2/pi) integral_4^6 rho(s)/(s-2)^3 ds.
The coefficient 2/pi is independently obtained by differentiating
the crossing-even kernel
2v^2/[pi(s-2)((s-2)^2-v^2)] at v=0 and dividing by two.

The exact positive weight integrals are
integral_4^6(s-2)^(-3) ds=3/32 and
integral_5^6(s-2)^(-3) ds=7/288.
Therefore

I_[4,6] > (1/2) 5lambda^2 (7/288)
          =35lambda^2/576 >lambda^2/20,
I_[4,6] < (2/3) 34lambda^2 (3/32)
          =17lambda^2/8 <3lambda^2.

These are bounds on the actual computed one-loop light cut, not on
an assumed constant spectral density. The exact angular expression
is retained in the definition; the inequalities justify its enclosure.

Let E be S6.113's complete one-loop b2 error majorant. Then

b2_tree+one_loop-I_[4,6] >4lambda-E-3lambda^2>0.

The last comparison is native exact rational arithmetic, and even
E+3lambda^2<10^-6(4lambda) still holds. No loop correction or cut
contribution is discarded merely because it is small.

This computes a positive finite-order improved functional. Calling it
a full improved positivity verdict would additionally require the
uncomputed higher-loop remainder and global dispersion hypotheses.
The finite integral alone cannot establish either.
