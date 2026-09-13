# Existing reference-profile split and both nonlinear clock contacts

Write Delta rho and Delta P for the full exact T-minus-S reference energy and pressure. Their zero-jet bounds are epsilon=10^6 B n^-7/2 and their required first five derivatives are integrable. By the exact KG Ward identity,

(Delta rho)'+3H(Delta rho+Delta P)=0.

On the reference slab, the existing S240 profile can be split algebraically into its comparison and selection contributions. The latter uses the SAME full localizer

T(X)=X^1024/[X^1024+(1-X)^1024]

times[-Delta P-(Delta rho+Delta P)(X-1)/2], with the common normalization by kappa0. Literal derivatives check T(1)=1 and T'(1)=T''(1)=0. This is not a new candidate, physical state or finite counterterm; the prescribed sum remains unchanged. No global smooth auxiliary profile is inferred.

## Entire physical Hessian and matching

For ADM N=1+n_lapse and h=a² exp(Q), tau=tr Q, the exact profile mixed density before normalization is

a³[-(Delta rho+Delta P)n_D n_G
+(Delta rho/2)(n_D tau_G+n_G tau_D)
-(Delta P/4)tau_D tau_G].

The full reference Gaussian state-difference first current is
a³[-Delta rho n_D+(Delta P/2)tau_D].
It is independently obtained from the isotropic covariance and ALL first Hamiltonian features. The profile's first current cancels it exactly; the shift mean vanishes by momentum reflection. The Ward identity also supplies the clock equation.

For |Delta rho|,|Delta P|<=epsilon, |tau|<=v and |n|<=v, the entire mixed coefficient is at most(2+1+1/4)a³epsilon<10^3epsilon. This gives the full profile bound10^9 B n^-7/2 v(D)v(G).

## Nonlinear common-clock map, not profile-only cancellation

The exact physical-to-affine metric relation is
zeta=v-(1/4)log R(t,N^-2).
Its first lapse derivative is delta and its second is4delta²-3delta, with delta=1/[2(1+t²)^3]. Therefore Q's second scalar chart derivative has trace6(4delta²-3delta)n_D n_G.

The state-difference current has the NONZERO second-chart contact

+3a³ Delta P(4delta²-3delta)n_D n_G.

The fixed profile has the opposite nonzero contact. Its whole same-clock mixed density is

a³[2DeltaJ n_D n_G+3Tc(n_D v_G+n_G v_D)-9Delta P v_D v_G],

DeltaJ=(21delta²-3delta)(-Delta P)/2
+(1-6delta)(-(Delta rho+Delta P)/2),
Tc=Delta rho-3delta Delta P.

A literal complete-profile lapse/clock expansion checks this density and the second-variation chain rule. The contacts cancel ONLY in the entire unprojected matched summand. At delta=1/2 the profile contact is positive3a³Delta P n_D n_G/2 and is generally nonzero.

## Required finite-cutoff remainder

At finite K, retain the original FULL profile but let Delta P_K be the actual projected current's pressure. The total nonlinear chart contact is then

-3a³(Delta P-Delta P_K)(4delta²-3delta)n_D n_G.

This is explicitly retained in the cutoff approximation. The mean tail has the full bound10^6 B K^-7, |4delta²-3delta|<=1, and a³<=64, so its contact is below10^9 B K^-7. It fits the same-clock regulator-tail allowance. Declaring the projected mean already matched would remove a genuine term and is rejected.
