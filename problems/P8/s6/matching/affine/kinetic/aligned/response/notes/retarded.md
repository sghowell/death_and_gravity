# A prepared retarded estimate for the actual leading vector block

The physical action is S6.42, unchanged. The calculation is at second
order in light perturbations and zero heavy initial data. Because the
stationary shifted source starts at second order and the vector starts
at second order in this state, coefficient variations multiplying the
vector enter its equation only at third and higher orders. Thus the
second-order forced operator is exactly the S6.42 quadratic operator.
This order-counting statement is not a bound on the higher orders.

## Derive the force before estimating it

Let S denote the second-order temporal source and sigma the longitudinal
vector potential, q>0. Retain the vector temporal variable t:

    L/a³=(t-S)²/2-q*sigma²/2+r*(sigma'-t)²/2, r=zeta*q.

The full temporal Euler equation gives t=(S+r*sigma')/(1+r), hence

    L/a³=r*(sigma'-S)²/[2(1+r)]-q*sigma²/2,
    (g²*sigma')'+a³*q*sigma=(g²*S)',
    g²=a³*zeta*q/(1+zeta*q).

Set v=g*sigma. The actual canonical equation is

    v''+W*v=J,
    W=q+1/zeta-U, U=g''/g,
    J=(g²*S)'/g=g*S'+2g'*S.

The source cannot be replaced by S or by g*S' alone. All primes hold
comoving momentum fixed. At this perturbative order S has no spatial
component, so neither transverse polarization is forced. At k_com=0
the temporal equation is t=S and all three spatial components remain
zero with zero data; no longitudinal normalization is inverted there.

## Frequency derivatives and the retarded energy bound

For 0<k_com²<=1, q<=1 at all times since a>=1. The S6.42 frequency
proof gives |U|<=15 and W>=D=1/zeta-15. Its alpha,beta expressions
give, with r>=0,

    r*alpha_r=r/(1+r)² <= 1/4,
    |r*(r*alpha_r)_r| <= 1/4,
    |r*beta_r|=|r*(5r-1)/(1+r)³| <= 1,
    |r*(r*beta_r)_r| <= 3.

The last two inequalities follow from the positive polynomial
majorants recorded in derivatives.py. With |H|<=2, |H'|<=4,
|H''|<=12 and |H'''|<=72, differentiation using r'=-2Hr yields

    |alpha'|<=1, |alpha''|<=6, |beta'|<=4, |beta''|<=56,
    |U'|<=74, |U''|<=688,
    |W'|<=78, |W''|<=712.

For completeness, H''=8u(u²-3)/(1+u²)^3 is bounded using
|u|/(1+u²)<=1/2 and (u²+3)/(1+u²)²<=3. The H''' numerator
has |u^4-6u²+1|<=u^4+6u²+1<=3(1+u²)². The alpha second
bound uses |1-r|<=1+r. These prove continuous bounds, not sampled ones.

For y''+W*y=F define E_y=sqrt(y'²+W*y²). Direct differentiation
gives (E_y²)'=2y'F+W'*y² and therefore

    E_y' <= |F|+39*E_y/D.

At zeros use sqrt(E_y²+epsilon²) and take epsilon down to zero.
On a forward interval of length <=1, with zero data, Gronwall yields
E_y<=exp(39/D)*||F||_infinity. For 0<zeta<=1/20000,
D>=19985>141² and exp(39/D)<=D/(D-39)<21/20. This is an
initial-value inverse estimate with specified data, not a spectral claim.

## Explicit local-response error

Suppose J is C², J=J'=0 initially, and |J|,|J'|,|J''|<=A. The
quasistatic Q=J/W has the same zero position and velocity data, and

    Q''=J''/W-2J'W'/W²-JW''/W²+2J(W')²/W³,
    |Q''|<=A*(1/D+868/D²+12168/D³).

The error e=v-Q solves e''+W*e=-Q'' with zero data. Thus

    |e| <= (21/20)*A/sqrt(D)*(1/D+868/D²+12168/D³)
         <= A/(100D).

The last inequality follows monotonically from D>=19985 and
sqrt(D)>141; the exact rational comparison is replayed. Since
|q-U|<=16, |J/W-zeta*J|<=16*zeta*A/D. Combining the two bounds,

    ||v-zeta*J||_infinity <= zeta*A/80.

This is an absolute norm error, not an error relative to J(u) at each
point. In particular, no division by zeros of the source occurs. A
bound on forcing amplitude alone is not the declared C² jet bound.
Unprepared J(0)=1 fails the proposed bound at the initial point, and
nonzero homogeneous heavy data are not source-proportional errors.

## Physical source and vector readouts

Restrict the interval to [-1/2,1/2]. Then 1<=a<=25/16. Let k=sqrt(k_com²).
If |S|,|S'|,|S''|,|S'''|<=E and S=S'=S''=0 initially, J=J'=0
and the above preparation applies. The actual normalization obeys

    g<=5*sqrt(zeta)*k/4, g>=sqrt(zeta)*k/sqrt(1+zeta),
    |g'/g|<=3, |g''/g|<=15, |g'''/g|<=74+15*3=119.

From J=g*S'+2g'*S,

    J'=g*S''+3g'*S'+2g''*S,
    J''=g*S'''+4g'*S''+5g''*S'+2g'''*S.

All three jet norms are bounded by A=(815/2)*sqrt(zeta)*k*E.
The physical spatial-vector amplitude is sqrt(q)*sigma. Since
sqrt(q)/g=sqrt((1+r)/(a³*zeta))<=101/(100*sqrt(zeta)), its local
approximation and error are

    W_spatial,local=zeta*sqrt(q)*(S'+2rho*S), rho=g'/g,
    ||W_spatial-W_spatial,local||_infinity <= 6*zeta*E.

For the temporal readout t-S=r*(sigma'-S)/(1+r), use the direct
energy estimate |v'|<=21A/20, |v|<=21A/(20sqrt(D)). Dividing by
g with the lower bound above gives

    ||t-S||_infinity <= 500*zeta*E.

The constants are conservative exact majorants. Both normalized
vector readouts acquire physical factor 1/tau. No arbitrary heavy
state, full nonlinear remainder, induced retarded quantum action,
loop error, interacting cutoff or V/G/B UV conclusion is included.
