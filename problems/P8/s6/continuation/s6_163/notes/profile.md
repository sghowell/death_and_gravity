# Global real bounds, local vacuum jets and a nonidentity control

Set x=z/R with R>0. The unique positive real eighth root
gives a real-analytic function on the whole real axis:

    f_R(z)=R x(1+x^8)^(-1/8),
    f'_R(z)=(1+x^8)^(-9/8),
    f''_R(z)=-9x^7/[R(1+x^8)^(17/8)].

Since |f_R|^8/R^8=x^8/(1+x^8)<1, the absolute argument
is strictly below R. The first derivative is positive
and at most one. For |x|<=1 the second derivative's
numerator has modulus at most nine; for |x|>=1 use
(1+x^8)^(17/8)>=|x|^17. Thus |f''_R|<=9/R.

Let Yupper be the unchanged rational upper bound on y^2.
At mF=10^200 and R=10^300,

    Yupper R^2/mF^2 <10^-4.

Consequently both masses are strictly in
(0.99mF,1.01mF) for every real z. This remains a pointwise
mass statement if z is replaced by any real scalar
expression; it is not a global Dirac propagator theorem.

The vacuum germ is

    f_R(z)=z-z^9/(8R^8)+9z^17/(128R^16)
             -51z^25/(1024R^24)+... .

The Taylor branch is analytic on |z|<R. Eight branch
points occur at z/R=exp[i(2j+1)pi/8]. Global real
analyticity is not complex entire analyticity.

The functions are not identical beyond the protected
jets. For the paired Euclidean log-determinant integrand

    log[(q^2+(m+y f_R)^2)(q^2+(m-y f_R)^2)],

its first difference from the original linear argument
has scalar-degree-ten coefficient

    -Y(q^2-m^2)/[2R^8(q^2+m^2)^2].

At q^2=0 this is Y/(2R^8 m^2), not zero. The expression
is a diagnostic of the paired logarithm before its Dirac
multiplicity and loop integration, not an extra term
silently added to a previously matched lower-point result.
