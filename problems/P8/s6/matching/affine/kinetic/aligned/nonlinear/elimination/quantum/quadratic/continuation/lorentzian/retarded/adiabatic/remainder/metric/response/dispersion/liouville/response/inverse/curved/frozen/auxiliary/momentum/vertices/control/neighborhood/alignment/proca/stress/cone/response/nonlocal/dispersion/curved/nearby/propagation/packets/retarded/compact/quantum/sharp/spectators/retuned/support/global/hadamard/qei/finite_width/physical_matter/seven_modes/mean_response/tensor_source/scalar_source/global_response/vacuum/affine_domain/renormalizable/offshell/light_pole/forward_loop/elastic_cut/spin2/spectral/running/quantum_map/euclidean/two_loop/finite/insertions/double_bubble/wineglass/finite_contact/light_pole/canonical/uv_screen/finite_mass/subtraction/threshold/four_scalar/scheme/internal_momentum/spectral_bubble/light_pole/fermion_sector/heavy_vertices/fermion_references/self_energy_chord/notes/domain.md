# Global complex gap and the momentum-dependent radius

Let q be any real Euclidean loop momentum and put
S=mF^2+q^2. Allow p=q+delta with
||delta||_1<=sqrt(S)/10. Then

    Re(p^2)>=q^2-21S/100,  |p^2|<=2S.

For b=0 or 1, mF>=2, the full parameter denominator
Delta=x mF^2+(1-x)b+x(1-x)p^2 obeys

    Re Delta >=(79/100)x[mF^2+(1-x)q^2]+(1-x)b.

The difference from the bound obtained by inserting
Re(p^2)>=q^2-21S/100 is the positive term
21x^2 mF^2/100. Thus for x>0, with t=q^2/mF^2,

    x/2 <= |Delta/mF^2| <=2(1+t),
    Re Delta >0.

The upper bound follows from b<=mF^2 and
x(1-x)<=1/4; its polynomial remainder is checked.
The principal logarithm has argument magnitude
below pi/2<2. Since log(2)<1,

    |log(Delta/mF^2)|<=log(1+t)-log(x)+4.

The x endpoint is controlled by
integral_0^1 -log(x) dx=1. The integrated logarithm
is bounded by log(1+t)+5; we use the looser +6.
This also controls the massless gauge denominator:
the endpoint singularity is integrable and its
momentum-dependent part is holomorphic.

The literal Euclidean Clifford matrices give
D(q)^dagger D(q)=S identity, hence
||S_F(q)||=1/sqrt(S). Moreover
||slash(delta)||<=||delta||_1. The Neumann identity
therefore gives ||S_F(q+delta)||<2/sqrt(S).
For the self-energy numerator
||slash(p)||<=|q|+||delta||_1<2sqrt(S), not a
Hermitian assumption about complex momenta.

The scalar kernel is consequently bounded by
3Y sqrt(S)[log(1+t)+6]/Q. The gauge logarithmic
and affine terms are bounded by
a Cf sqrt(S){8[log(1+t)+5]+4}/Q, which is at most
12a Cf sqrt(S)[log(1+t)+6]/Q. Together,

    ||Sigma_MS(p)|| <= (3Y+12a Cf)/Q
                       sqrt(S)[log(1+t)+6].

For the forward unit disc |s-2|<=1, u=4-s,t=0,
choose analytic center-of-mass momenta with
E^2=s/4 and P^2=s/4-1. Each Euclidean incoming
momentum has L1 norm below two. All masses are
one and their sum is zero. Every loop-route
prefix contains at most three such momenta,
so its L1 norm is below six, hence below 18.

Scale all external momenta by a complex variable
zeta. For each q choose R(q)=sqrt(S)/360.
On |zeta|<=R(q) every route shift is below
sqrt(S)/20, strictly inside the proved domain.
For mF>=720, R(q)>=2 uniformly. No expansion
of q/mF or upper momentum cutoff is introduced.

The positive inner parameter gap and the outer
Dirac resolvent bound establish the continuation
for all real q. No unproved translation of a
real integration contour at complex momenta is
used. This two-point-kernel argument does not
automatically apply to a vertex subgraph with
an independent scalar-transfer threshold.
