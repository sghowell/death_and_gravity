# Regular lower dictionary and the unavoidable boundary change

In physical X variables, the actual zero-braiding dictionary
is the regular-singular coefficient ODE

    q_X+[1/(2X)-3R_X/(4R)]q=3R_X R_u/(4R)=g(u,X).

For X>0 its integrating factor is I=sqrt(X)R^-3/4.
The solution regular at zero is most usefully written

    q(u,X)=X R(u,X)^(3/4)
       int_0^1 sqrt(s) R(u,sX)^(-3/4) g(u,sX) ds.

This formula is real on both sides of zero. The segment
sX stays inside the connected positive-R domain. On every
compact parameter subdomain, all integrand derivatives
are bounded by an integrable multiple of sqrt(s), and
the positive fractional-power branch has a common analytic
neighborhood. It therefore defines a real-analytic q.
This is an integral over coefficient space, not over physical
time, a retarded kernel or a quantum state history.

Since R_X=-2nX/h+O(X^2) and R_u=n h_u X^2/h^2+O(X^3),

    g=-(3n^2 h_u/(2h^3))X^3+O(X^4),
    q=-(n^2 h_u/(3h^3))X^4+O(X^5).

The homogeneous alternative C(u)R^3/4/sqrt(X) cannot be
regular across zero unless C vanishes. This proves uniqueness
among analytic vacuum-regular lower dictionaries.

The actual source-convention coefficients are

    J3=(q+R_u/2)/(4pX),
    J2=F-X q_u+3X(q+R_u)^2/(16p^2).

The literal affine source equations give Q1=q, Q2=-2q_X
and P+X q_u=F. The covariant divergence of q d(phi)
then removes the braiding exactly. Every lower term and
the sign of the scalar boundary contribution are retained.
At zero, J3/X -> n h_u/(4h^2); J2 is analytic as well.

The old clock condition q(u,1)=0 cannot be retained.
For R=1+B(X)(X-1)/h define

    W(R)=(4/3)R^-3/4+4R^1/4-16/3,
    W'(R)=(R-1)R^-7/4.

Writing y=R^1/4 gives the exact factorization

    W=4(y-1)^2(3y^2+2y+1)/(3y^3)>=0,

with equality only at R=1. Since R_u=(h_u/h)(1-R),
I g=-(3h_u/(4h))sqrt(X) d_X W(R).
Integration by parts at R(0)=R(1)=1 gives

    q_reg(u,1)=(3h_u/(8h)) int_0^1 W(R(u,z))/sqrt(z) dz.

For finite u and 0<z<1 the integrand is positive. Thus this
boundary has the sign of u and is nonzero for u!=0.
The old boundary instead leaves its negative as the
coefficient C(u) of the divergent homogeneous term. Although
C(0)=0, C'(0)=-(9/4)int_0^1 W(R(0,z))/sqrt(z) dz<0.
The isolated u=0 cancellation does not provide a regular
joint vacuum neighborhood. The smooth plateau version
has the same argument with vanishing integrand below 1/4.

There is also an explicit global bound. On 0<X<1, set
X=(1+tanh(y))/2. Then T'=n cosh^2(y)/cosh^2(ny)<=n;
outside that interval on the stated domain, the elementary
power ratio also bounds |T'| by n. Positive exponential
Taylor terms give |w'|<=4sqrt(n)<=4n. Consequently
|R_X|<7n/h and |R_u|<5|h_u|/(4h^2).
Using 1/2<R<=6/5 and int_0^1 sqrt(s) ds=2/3 yields

    |q(u,X)| <32n |h_u|/h^3 <=96n/h^2.

This bound covers all real u, including both tails.
It does not supply derivative/loop or heavy-threshold
bounds for a propagating completion.
