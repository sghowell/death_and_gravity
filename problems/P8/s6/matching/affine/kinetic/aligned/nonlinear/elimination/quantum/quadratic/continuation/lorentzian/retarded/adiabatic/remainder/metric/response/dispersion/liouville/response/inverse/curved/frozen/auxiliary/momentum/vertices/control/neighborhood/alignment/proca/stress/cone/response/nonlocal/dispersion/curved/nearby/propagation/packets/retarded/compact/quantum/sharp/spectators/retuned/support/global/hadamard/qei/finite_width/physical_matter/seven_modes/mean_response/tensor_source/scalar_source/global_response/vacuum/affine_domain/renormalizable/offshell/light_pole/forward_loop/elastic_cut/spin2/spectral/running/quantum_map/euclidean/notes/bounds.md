# Actual rational bounds without floating high-mass cancellation

Let A=g/(288M), using the frozen actual polynomial parameters.
The same quantity already bounds the finite light kinetic
counterterm in S6.112. Exact rational checks give

    0<A<3 times 10^-208,
    A/(1-A)<3 times 10^-208,
    A(2-A)/(1-A)^2<6 times 10^-208.

The known-kernel tails beyond one retained insertion obey

    A^2/(1-A)<10^-414,
    A^2(3-2A)/(1-A)^2<10^-414.

These are dimensionless relative propagator factors, not
absolute scattering coefficients or all-higher-loop errors.

There are also nontrivial lower bounds at every finite y.
For x in [1/4,3/4], M>=1 implies d(x)<=M and
x(1-x)>=3/16, so b>=bmin=3/(16M).
Using the positive integral representation and pi<4 gives

    Q(y) > g bmin^2(y+1)/[1024(1+bmin(y+1))] > 0.

The global upper and the quadratic small-increment bound combine as

    Q(y) < min(A, g(y+1)/(864M^2)).

The lower estimate integrates only the middle half of the actual
parameter interval; it does not cut off loop momenta. Both bounds
are valid for arbitrary finite nonnegative rational inputs and
follow from continuous inequalities for all real y>=0.
