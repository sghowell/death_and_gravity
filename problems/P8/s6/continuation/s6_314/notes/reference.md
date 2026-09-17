# Unexpanded reference comparison and two-marked dressing

For a marked state of total energy R<=x, combine the new finite
conversion modulus with S300's index bound and the Gamma-product
derivative bound:
 |delta Delta|<5000R(1-lnR)/kappa,
 |delta a|<40R/kappa,
 |delta ln F(a)|<480R/kappa^2,
 F(a)=exp(-gamma_E*a)/Gamma(1+a).
Since |ln x|<=|ln R|,
 |ln[P_sigma(x)/P0(x)]|
 <6000R(1-lnR)/kappa<1/2.
Here R(1-lnR)<=4/8 because R<=1/8 and ln8<3.
It follows that
 |P_sigma(x)/P0(x)-1|<12000R(1-lnR)/kappa.

Put z=1-R/x. Both physical indices are nonnegative, so the real
mean-value theorem gives
 |z^a_sigma-z^a0|<=|delta a|*|lnz|.
Writing the remaining-energy ratio as
 [P_sigma(x)/P0(x)]z^a_sigma-z^a0
then yields
 |P_sigma(x-R)-P0(x-R)|/P0(x)
 <[12000R(1-lnR)+40R|ln(1-R/x)|]/kappa.
The angularly integrated elastic seed is a0*dR/R,
a0<4/(5kappa). Using the exact integrals
 int_0^x (1-lnR)dR=x(2-lnx),
 int_0^x -ln(1-R/x)dR=x,
we obtain the convenient conservative bound
 |E1|/P0<10^5*x*(2-lnx)/kappa^2.
It is uniform down to arbitrarily small x and tends to zero.
The resolution power has never been expanded in a0*lnx.

## Two-marked finite signed seed

The full r2 measure is finite by notes/matching.md. S309's additional
leading-soft signed-series proof applies on its fixed two-marked state
sigma_ab because S300/S301 bounds depend on TOTAL marked energy.
For n=N+2 labels, summing the unordered marked pairs on the common
energy simplex gives
 choose(N+2,2)/(N+2)!=1/(2!*N!).
The seed measure already includes1/2!, so the extra emissions have
only their1/N!. Their total remaining energy is x-a-b, not x.
No two independent individual-energy cuts are substituted.

The common absolute series bound by6 permits dominated convergence
against the finite signed seed at fixed positive x. The physical
kernel is nonnegative, and
 P_sigma_ab(x-a-b)/P0(x)
 =[P_sigma_ab(x)/P0(x)]*(1-(a+b)/x)^a_sigma_ab
 <=1+4250/kappa<2.
Therefore
 |D2|/P0<2*TV(R2;x)<4*10^-725*x+4*10^-652*x^2.

## The named reference and its precise strength

Define P_match2=P0+D1+E1+D2, with D1 exactly the unchanged S309
one-marked remainder dressing. Triangle bounds give the expression
in the formulation. At x<=1/8, its three contributions are bounded by
 2*10^-768,
 10^5*(5/8)/kappa^2,
 4*10^-725/8+4*10^-652/64.
Their exact rational sum is<10^-653. Each has a separate uniformly
vanishing majorant at x->0. Since P0>0, P_match2 is pointwise positive.

This statement is not event-measure positivity, monotonicity in x,
unitarity or full inclusive matching. The two-marked r2 is dressed
on its defined full sigma_ab; its internal counterterms are not
promoted to actual higher-real matrix elements. Higher nonleading
seed sectors and their compatibility still need derivation.
