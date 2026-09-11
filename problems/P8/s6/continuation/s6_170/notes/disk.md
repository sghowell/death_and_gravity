# Variable complex disks and exact derivative control

Put L=sqrt(t^2+tau^2), D=1+t^2, q=p/D^2 and E=sqrt(q^2+m0^2).
Assume 0<tau<=1 and 0<=Delta/m<=.003. On |z-t|<=L/1024,

    |(1+z^2)/D-1|<=2/1024+1/1024^2<1/256.

This uses L^2<=D and |t|L<=D. Thus the complex momentum q(z)
is close to q and |q(z)^2-q^2|<q^2/16.

For x=t/tau the same disk has radius sqrt(1+x^2)/1024,
at most max(1,|x|)/512. With r=1/512 the binomial bound gives

    |z_x^8-x^8|/(1+x^8)<=8r(1+r)^7<1/32.

The positive-real eighth-root branch extends to this disk.
Consequently |s(z_x)|<2 and |s'(z_x)|<2w(x),
where w(x)=(1+x^8)^(-9/8) and integral_R w=2.
The complex mass differs from m by less than2Delta. Therefore

    |M(z)^2-m0^2|/m0^2
       <=[4(.003)+4(.003)^2+1-.99^2]/.99^2 <1/24.

Combining both errors gives
|omega(z)^2-E^2|/E^2<1/16+1/24<1/9. The branch is nonzero,
analytic, and |omega(z)|>.9E.

Also |H(z)|<5L/D, |q(z)|<2q and |M(z)|<2m. These inequalities
bound the exact initial connection by

    G0=16q/E^2 [Delta/tau w(t/tau)+mL/D].

Since sqrt(1+x^2)w(x)<=2, L^2/D<=1 and q/E<=1,
G0/E<=32/(LE).

Fix any positive integer N. Set K_N=8192N and nested disk radii

    rho_j=1/1024-j/(2048N), j=0,...,N.

Suppose LE>=2K_N. On disk j the induction is

    |e_j|>=c_j E, c_j=.9(1-1/(4N))^j>.5,
    |g_j|<=G0 (K_N/(LE))^j.

The ratio |g_j/e_j|<=64/(LE)<=1/(256N) is small.
The exact square root sqrt(1+(g_j/e_j)^2) has modulus at least
1-1/(4N). Bernoulli's inequality gives
c_j>=.9(1-j/(4N))>=.675, for every j<=N.
For |v|<=1/2, the analytic atan branch obeys |atan v|<=2|v|.
Cauchy's derivative estimate with radius loss L/(2048N) then yields

    |g_(j+1)| <=4096N G_j/(LE)<=G_j K_N/(LE).

All estimates are taken on disks with the same fixed real center;
no unproved comparison of weights at different centers is needed.
The lower square-root and atan branches agree with the real exact
positive-frequency rotations. The final radius is still1/2048.

For twenty frames the actual tau*m0>2K_20 holds at every time and
momentum. It is NOT asserted for unbounded N. Instead, for |t|<=1,

    LE>=tau*p/4,

and for |t|>=1,

    LE>=max(m0|t|,p/(4|t|^3))>=(m0^3p/4)^(1/4).

Thus any fixed N is allowed at all times once
p>=max(8K_N/tau,4(2K_N)^4/m0^3). This order-dependent high-momentum
threshold is what supports the later all-order regularity argument.
