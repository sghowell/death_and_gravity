# Uniform canonical vertex bounds and the remaining propagator gap

These bounds are for coefficients, not inclusive probabilities.
Use complex Euclidean Frobenius norms for fields and connections,
the induced Euclidean operator norm for inverse metric tensors,
and Euclidean norms for momenta. Multiplication by eta is an
isometry for these norms.

For a full mixed log-density coefficient,
 |tr(eta A1...eta Ar)|<=4 product||A_i||F.
The full permutation sum therefore gives
 |L_r|<=2^(r+1)(r-1)! product||A_i||F.
Its scalar exponential generating majorant is
 sum_(r>=1)2^(r+1)z^r/r=-2 ln(1-2z).
The anchored exponential recursion consequently bounds the density
by(1-2z)^-2. The inverse permutation sum is bounded by(1-2z)^-1.
Subset convolution bounds the density inverse by(1-2z)^-3.
Their coefficients prove
 |D_r|<=r!2^r(r+1), ||I_r||op<=r!2^r,
 ||K_r||op<=r!2^r binom(r+2,2),
with the field-norm product restored in each expression.
These are analytic majorant identities, valid at every order,
not estimates inferred from the first few tensor evaluations.

The undoubled linear connection satisfies
 ||C1(A,p)||F<=3||p||2||A||F/2.
Each of its three terms is a permutation of p tensor A; eta has norm1.
In Gamma_r, sum over which of the r fields carries the derivative.
If every momentum norm is<=L this gives
 ||Gamma_r||F<=3*2^(r-1)r!*L product||A_i||F,
with generating majorant3Lz/(1-2z).

For B(K,G,H), the first trace contraction is bounded by
4||K||op||G||F||H||F: the4x4 Frobenius norm of K is at most twice
its operator norm and the single connection trace costs at most2.
The second contraction is at most2||K||op||G||F||H||F by matrix
Cauchy-Schwarz. Hence the conservative common constant8 is valid.
Multiplying the three generating majorants and the action factor1/2
gives36L^2 z^2/(1-2z)^5. Its rth mixed coefficient is
36r!2^(r-2)binom(r+2,4)L^2 product||A_i||F.
For r>=3, binom(r+2,4)<=2^(r+2), and36*4^r<=32^r;
this proves the simpler r!32^r L^2 envelope.
Restore kappa^(1-r/2).

A scalar vertex follows directly by a quadratic-form bound:
 |V_phi2,h^r|<=r!2^r
 [binom(r+2,2)||p||2||q||2+(r+1)m2]
 product||A_i||F/kappa^(r/2).
The heavy and contact potential vertices have their coupling magnitude
times the determinant bound. Internal field norms are not bounded
by1 in this statement.

## Arbitrary finite-N hard cuts versus pure-soft cuts

The S312 hard-cut geometry was written for two real gravitons, but its
proof uses only Qa,Qb future causal and Qa+Qb=Q, with total energy W.
Those conditions hold for every subset of any finite collection of
future null emissions. With A=p3+Qa,B=p4+Qb it gives unchanged
 -K^2>(tau+W^2)/300000
for a mixed hard Einstein cut at positive associated Born transfer.
Its Kallen and recoil inequalities depend on W, not on the number of
rays. The timelike hard range45/8<=K^2<=16, heavy inverse bound2/n,
and light-scalar gap3*w_subset/8 likewise use only future subset sums.
Zero-energy leaf cases are handled by external seeds, not internal
propagators.

This does NOT cover pure-soft graviton cuts Q_S^2. Three or more
soft rays can have nested overlapping collinear clusters. Their
complete vertex/contact cancellations need a quantitative conserved-
current bound uniform in those regions and with controlled growth
in N. Neither the finite pair result nor the graph/vertex majorants
here supply it. A small gravitational coupling alone is not a proof
that an unspecified all-N remainder is small.
