# A uniform analytic disk and four derivative-control steps

Fix real x and consider |z-x|<=r=1/64. The difference of
eighth powers and a split into |x|<=1 and |x|>=1 give

    |z^8-x^8|/(1+x^8)<=8r(1+r)^7<1/4.

Thus 1+z^8 stays strictly in the right half plane.
The branch s(z)=z/(1+z^8)^(1/8) extending the real profile
is analytic throughout this strip. In the disk centered at x,

    |s(z)|<2,
    |s'(z)|<2w(x), w(x)=(1+x^8)^(-9/8).

For the first inequality use
(|x|+r)/(1+x^8)^(1/8)<=1+r and the denominator factor
(4/3)^(1/8). For the second use (4/3)^(9/8)<2.
The real weight has integral two.

Let E=sqrt(p^2+m0^2), m0=.99m. Since Delta/m0<1/100,
the complex mass satisfies |M(z)-m|<=2Delta and

    |M(z)^2-m^2|<=4mDelta+4Delta^2<.1m^2.

The real part of p^2+M(z)^2 is therefore positive.
Its principal square root omega is analytic, has its
physical positive real value, and obeys |omega|>=.9E.
For the initial mass-diagonal frame,

    e0=tau omega,
    |g0|=|pDelta s'(z)/(2omega^2)|
         <=2pDelta/E^2 w(x).

For j=0,...,4 set

    r_j=1/64-j/512,
    c_j=.9(15/16)^j>1/2,
    G_j=2pDelta/E^2 [1024/(tau E)]^j.

We prove on the nested disk of radius r_j that
|e_j|>=c_j tau E and |g_j|<=G_j w(x).
The initial case was just established.

Suppose these hold at j<4. Write q_j=g_j/e_j.
Using w<=1, p<=E and E>=m0,

    |q_j|<=2G_j/(tau E)
         <=4Delta/(tau m0^2)[1024/(tau m0)]^j<1/4,

because tau*m0>=2^20. Hence the branch
e_(j+1)=e_j sqrt(1+q_j^2) is analytic and its square-root
factor has modulus at least sqrt(15/16)>=15/16.
This proves the next frequency floor.

On the smaller disk Cauchy's derivative estimate,
with radius loss1/512, gives

    |q_j'|<=512 G_j w(x)/(c_j tau E).

The exact next connection has absolute value
|q_j'|/[2|1+q_j^2|]. Its last denominator factor is
strictly greater than one. With c_j>1/2 this gives

    |g_(j+1)|<=1024 G_j w(x)/(tau E)=G_(j+1)w(x).

All four steps have positive disk radii and strict
complex gaps. The argument is not just a real-time
eigenvalue check and is not an expansion in Delta.
Every constant is explicit and independent of momentum
and real time under the stated parameter conditions.
