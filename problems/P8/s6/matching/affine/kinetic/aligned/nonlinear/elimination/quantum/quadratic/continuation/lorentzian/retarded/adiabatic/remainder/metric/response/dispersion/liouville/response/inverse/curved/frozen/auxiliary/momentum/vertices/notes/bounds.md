# Continuous time and momentum majorants

These are raw canonical phase-space estimates, prior to
free-mode normalization and the factors (M*tau)^(2-n).

## Uniform actual coefficient bounds

Each actual lapse derivative is rational in real u.
Its reduced denominator is an even polynomial with
nonnegative coefficients and strictly positive constant.
If c(u)=P(u)/D(u), |u|<=1/2, then

|c(u)| <= sum_j |P_j|*2^-j / D_0.

The verifier checks that denominator certificate for
every one of the 70 actual derivatives. It also retains
the independent S6.75 Cauchy bound

|partial_N^k c_alpha(1,u)| <= 10000*k!*100^k,

from its N-radius 1/100, unit nine-invariant polydisc
and complete Hamiltonian bound 10000. The smaller of
these two justified bounds may be used.

The common bounds for N orders zero through four are

736125197/13107200,
327/32,
3817522971573/53687091200,
3838025560951302003/2147483648000000,
256462152065474570379/6871947673600000.

They are bounds on monomial coefficients after combining
the actual terms, not on selected summands that might
have been dropped. The original primitive, Q source,
inverse Maxwell normalization and all 14 monomials remain.

## Positive degree algebra

Use a nonnegative series truncated after degree four.
At each degree the norm sums absolute coefficients over
all labelled masks. Ordinary positive convolution is an
upper bound on nilpotent multiplication because it also
counts products with repeated labels, which actually vanish.

If every individual phase-component amplitude of each
leg is bounded by S and there are at most four legs,
f=4S*t is a common linear seed. Let K bound the modulus
of every coordinate derivative on the Fourier masks used.
Let D bound 1/|k| on the nonempty proper masks only.

A componentwise metric perturbation bound is 3f.
Neumann inversion gives

g<=1+3f,
g^-1<=1+sum_{r=1}^4 3^(r-1)*(3f)^r.

The six determinant products give a positive bound on
det(g)-1 after deleting its exact zero constant.
Absolute finite binomial series then bound sqrt(g) and
its inverse. Christoffel and scalar-curvature counts are

Gamma <= (9/2) K*g^-1*(3f),
R <= 9*g^-1*(6K*Gamma+18*Gamma^2).

On I, |H|<=2 and |ell|<=1/10. Thus every free momentum
entry is bounded by 2+6f. At each York degree the source
is bounded by

3K*pi+9Gamma*pi
 +(3/2)g^-1*[K*(1/10+f)*f+9K*f^2].

The last term retains both vector pieces:
Pi^j F_ij contributes at most 6K*f^2,
and -W_i div(Pi) contributes at most 3K*f^2.

The exact full tensor inverse bound of S6.76 implies

|(L E^-1 source)_ij| <= sqrt(6)D*max_i|source_i|
 <3D*max_i|source_i|.

Add that degree's bound to pi, and repeat through cubic
order. This is a uniform nonzero-transfer estimate;
no denominator at the zero total momentum is introduced.

For mixed=pi*g/sqrt(g), an entry bound is
3*pi*g/sqrt(g). The nine physical invariant majorants
are obtained with their actual index counts. In
particular the electric term carries 1/zeta, magnetic
term carries zeta, and the full vector divergence has
three derivatives. Only exactly absent background and
linear coefficients of squared invariants are removed.

## Stationary and moving-boundary counts

Let B_k be the five common coefficient bounds above.
A bound for L_k is B_k times the sum of all nine
invariant majorants. A bound for Q_k is B_k times

dp^2+dp*j+j^2+dc^2.

The 14 actual monomials contain exactly these quadratic
powers. Their coefficients already include cross-term
multiplicities.

Use |1/A2|<20, hence |1/(2A2)|<10.
Positive counterparts of n1,F2 and every stationary H
term in notes/reduction.md give the raw density bound.
Multiply by the complete volume majorant.

The moving term 2H*pi:g contributes at most 36*pi*g.
Also |H'+3H^2|<=16, so its explicit boundary contributes
at most 96f+120f^2; the matter term adds f/10.
Keeping these positive bounds is safe even where
spatial TT integration makes their higher coefficients
vanish exactly.

This proves the H3/H4 coefficient bounds for all fields,
all polarizations and every momentum configuration
satisfying the stated amplitude/derivative/transfer
bounds. It is not established by the finite fixtures
alone. The mixed exact fixtures are independently
checked to lie below the S=1,K=4,D=1 instance.

## Why this is not yet a physical cutoff

S bounds raw phase components, not normalized oscillator
columns. The latter require the actual coupled scalar
kinetic matrices and moving momentum boundary, both
tensor polarizations and transverse/longitudinal Proca
normalizations. Finite-time integration, exchange,
species multiplicities and momentum-fiber measures
are also absent from the raw coefficient bound.

The raw estimate supplies a quantitative input to those
next steps. It is neither a transition probability nor
a Wilsonian UV or all-orders remainder bound.
