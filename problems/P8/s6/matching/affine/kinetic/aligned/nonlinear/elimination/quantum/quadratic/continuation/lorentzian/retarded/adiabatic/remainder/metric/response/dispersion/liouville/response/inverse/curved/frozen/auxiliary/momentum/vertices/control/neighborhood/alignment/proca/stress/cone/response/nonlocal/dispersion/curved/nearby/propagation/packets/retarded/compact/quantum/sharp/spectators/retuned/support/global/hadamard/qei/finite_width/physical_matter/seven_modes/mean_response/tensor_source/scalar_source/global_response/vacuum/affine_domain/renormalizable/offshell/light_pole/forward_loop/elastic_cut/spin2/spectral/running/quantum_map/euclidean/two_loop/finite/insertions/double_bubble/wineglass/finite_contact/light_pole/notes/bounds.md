# Actual two-loop bound and pole/residue conclusion

If F is holomorphic on a neighborhood of |s-1|<=R and |F|<=S,
its Taylor coefficients obey |a_n|<=S/R^n. Subtracting a_0
and a_1(s-1), then summing the tail, gives

    |OS F(s)| <= S |s-1|^2/[R(R-|s-1|)].

With R=2 and |s-1|<=1 the coefficient is at most S/2.
The mixed integral uses its already finite constant-reference
difference; OS annihilates that reference. The same-heavy,
finite-sunset and nested-R integrals are bounded before OS.
The local sunset instead uses its directly integrated second
derivative and the line Taylor factor one half.

All positive absolute two-loop coefficients are therefore

    local sunset:       L^2/(3 Qnorm^2),
    mixed sunset:       (33/2)Lg/Qnorm^2,
    same-heavy sunset:  (5/4)g^2/Qnorm^2,
    finite sunset:      28g^2/Qnorm^2,
    nested R:           2g^2/Qnorm^2,
    nested alpha:       alpha_upper B1.

The local one-loop parameter/finite-contact pieces and the
new outer mass/kinetic terms have already been handled by OS;
they do not supply additional nonlocal two-point contributions.
In particular the alpha term keeps its subtraction-generated
overall divergence grouped with the outer on-shell operation.

At the exact inherited parameters,

    M=2lambda/gamma+2,
    g=1/67108864,
    B1=g/[864 M^2(1-1/M)],
    alpha_upper=g/(288M).

L is the actual positive polynomial quartic from S6.110,
not a newly selected coupling. Since pi>3, Qnorm^2>20736.
The exact rational bound is

    B2 = alpha_upper B1
       + [L^2/3+(33/2)Lg+(125/4)g^2]/20736
       = approximately 3.3463030014984704753 times 10^-19.

Every group upper bound is strictly positive and native exact
rational comparisons establish B1+B2<10^-18. No floating-point
comparison is used as a proof gate.

Define Pi_j,R(s)=(s-1)^2 H_j(s). Its two fixed OS zeros
make H_j holomorphic throughout the disc. The derived bounds
give |H1+H2|<=B1+B2. For the inverse propagator through
two loops,

    D_[2](s)=(s-1)[1+(s-1)(H1(s)+H2(s))].

On the unit disc the bracket differs from one by less than
10^-18, hence is never zero. The sole zero of D_[2] is
s=1, it is simple, and its derivative is one. The propagator
has unit residue there. The same reasoning holds for every
exact rational inner radius 0<r<=1; the factored error is
bounded by r(B1+B2). Unsupported/inexact radii are rejected
before any cache lookup.

This is a statement about D_[2], not the exact all-orders
propagator. Dropping heavy-mass suppression where convergence
survives makes B2 much larger than the very tight inherited B1.
Unequal conservative upper bounds do not establish perturbative
growth or the sign of the actual correction. No all-energy,
global-pole, Regge or nonperturbative conclusion follows.
