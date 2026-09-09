# Explicit replacement of the fixed scalar profiles

S6.81 retained the old fixed scalar profile solely to compare
the changed action. Its old stress cancellation did not transfer.
This candidate explicitly replaces it. From the new action,
the specified complete clock and the SAME fixed state, compute
the new finite functions rho_new(u), p_new(u) once by the
subtracted continuum integrals and matched local terms in
integrands(). These are c-number coefficient functions.
Thereafter hold them fixed under metric, clock and state
variation; dependence on the clock argument is varied as for an
ordinary local coefficient, not recomputed as a stress functional.

On the clock tube set

    P_new(u,x)=-p_new(u)+(rho_new(u)+p_new(u))*(x+1)/2,
    x=g^(mu nu) partial_mu(phi) partial_nu(phi).

Multiply by the same globally smooth window
w(x)=1-chi(64(x+1)^2). It is identically one for
|x+1|<=1/8 and zero for |x+1|>=sqrt(2)/8.
The original |x+1|<=1/10 tube is flat; an open neighborhood
of x=0 and all local jets there are untouched. This does
not establish existence of a healthy Minkowski vacuum.

The pressure is exactly p_old. Hence the literal new-minus-old
coefficient is

    Delta P = w(x)*(rho_new-rho_old)*(x+1)/2.

At x=-1, the local stress of P_new has energy -rho_new and
pressure -p_new. Its scalar equation on the clock is
rho_new'+3H*(rho_new+p_new), which vanishes by the new
ordinary Proca Ward identity. Therefore the retained Gaussian
one-point clock equations cancel exactly after this declared
profile replacement. The state construction is not circular:
P_new contains no vector, leaves its clock operator and initial
state unchanged, and is frozen before checking the mean equation.
No other field loops have been included.

For eta_j bounding both new stress functions and d_j bounding
their energy change, on the clock tube

    |partial_u^j P_new|<=1.1 eta_j, |partial_u^j P_new,x|<=eta_j,
    |partial_u^j Delta P|<=d_j/20,
    |partial_u^j Delta P_x|<=d_j/2,

and both second x derivatives vanish. On the global support
|x+1|<1/4, if W_k bounds the k-th window derivative, the
mixed derivative bounds for j+k<=5 are

    new: [(5/4)W_k+k W_(k-1)] eta_j,
    change: [W_k/8+(k/2)W_(k-1)] d_j.

The k-1 term is absent for k=0. The 21 entries for each
triangle are all below 10^-770 at L=10^400, including the
transition. The largest new-profile upper bound is less than
1.415*10^-771 and the change bound less than 2.723*10^-772.
These are continuum rational estimates, not sampled maxima.

For the actual point-chart density, x=-N^-2 and physical
volume is N*exp(3omega+3v), with
omega=-log((h-1+N^-2)/h)/4. The replacement has zero clock
value, lapse force Delta rho, zero scale force, and quadratic
coefficients

    n^2: (3/h-1)*Delta rho/2,
    nv: 3*Delta rho,
    v^2: 0.

All are derived directly. For h>=1 their absolute upper
bounds are d0,3d0,0. The full new profile's corresponding
bounds are 15eta0/8,15eta0/2,9eta0/2, from the same universal
point-chart coefficient identities as the parent. They too
are below 10^-770. These local contact coefficients do not
bound the nonlocal quantum light response.
