# Full low/high decomposition and positive moment Gram bound

Assume the explicit physical and differentiability premises of FORMULATION. Set w=S-2. The fixed-t dispersion kernel is

2v²/[pi (S+t/2-2)((S+t/2-2)²-v²)].

Taking its v² coefficient, its first transfer derivative and its v4 coefficient gives the three exact full relations in FORMULATION. The definitions use half the second and one twenty-fourth the fourth derivative, respectively.

Let B_low,B_high be the b20 contributions on[4,K] and(K,infinity), and J_low,J_high the corresponding inverse-fourth moments. Write

Tpos=(2/pi) integral rho_t/w³>=0.

Then b20=B_low+B_high and b21=Tpos-3(J_low+J_high)/2. Since w>K-2 on the high interval, J_high<=B_high/(K-2). Thus the following EXACT decomposition is nonnegative:

b21+3b20/[2(K-2)]+3J_low/2
 =Tpos+3B_low/[2(K-2)]
       +(3/2)[B_high/(K-2)-J_high]>=0.

An atom at K belongs to the low interval, so there is no double counting. The proof works for positive measures and retained heavy atoms, not only smooth spectral densities. The light cut below K is present throughout; K is an arbitrary analysis split.

Using the physical coefficient upper bounds
b20<=4lambda(1+delta0) and b21<=-3gamma(1-delta1)
gives
J_low>=2gamma(1-delta1)-4lambda(1+delta0)/(K-2).

At K10^198 and the named tolerances, this is at least gamma-8lambda/(K-2), strictly greater than gamma/5 by an exact rational comparison. No small full-loop error is used.

If a FUTURE calculation instead provides a full budget J_low<=epsilon gamma, and2(1-delta1)-epsilon>0, the simultaneous conditions imply

K<=2+4lambda(1+delta0)/{gamma[2(1-delta1)-epsilon]}.

This constrains the analysis split at which those matching and FULL-budget premises can coexist. It is not a physical cutoff, and no such full budget is established here.

For the independent higher-coefficient consequence, introduce the finite positive measure

dnu=(2/pi)rho dS/w³.

Then b20=integral dnu, J_total=integral dnu/w and b40=integral dnu/w². Its Gram matrix for the functions1 and1/w is positive semidefinite: integrating any real square proves positivity, or Cauchy-Schwarz gives

J_total²<=b20*b40.

The full b20 is positive when nonzero absorptive weight is present, as is forced by negative b21 under the named hypotheses. Since Tpos>=0,
J_total>=-2b21/3 for b21<0. Consequently

b40>=4b21²/(9b20)
    >=gamma²(1-delta1)²/[lambda(1+delta0)].

For the named errors this is gamma²/(8lambda)>0; for exact b20,b21 matching it is gamma²/lambda. These are necessary bounds on the physical FULL forward v4 coefficient, whereas the original tree coefficient is zero. They neither compute a loop counterterm nor permit an arbitrary physical coefficient assignment.

A positive two-atom fixture verifies the Gram determinant as p0*p1*(x0-x1)². A one-atom measure saturates it, as described in notes/positive-pole.md. Saturation in the cone of positive dispersion measures is not proof of attainability by an exact unitary interacting S matrix.
