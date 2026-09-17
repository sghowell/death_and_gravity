# Exact tensor majorants and their limited domains

## Two-ray field bound

By spatial rotation put q1 along z and q2 in the xz plane. Let a,b>0
be their energies and r the stereographic half-angle parameter.
Use the unnormalized plus/cross TT basis at each ray; each basis
tensor has Frobenius norm sqrt(2).

S311 computed all24 coefficients of the pair current contracted with
an arbitrary conserved upper root. Here each group of six coefficients
is compared to the spatial temporal-pair tensor entries

(Txx,2Txy,2Txz,Tyy,2Tyz,Tzz).

All24 identities hold exactly as rational functions of a,b,r.
The offdiagonal factors are retained: the sum of the absolute values
of these six coefficients is the entrywise l1 norm of the symmetric
spatial matrix, which bounds its Frobenius norm.

S311's coefficient budgets for the four plus/cross pairs are
530,128,88,268 times(a+b)^2/(a*b) at unit canonical coupling.
For a unit-Frobenius TT tensor A=alpha*plus+beta*cross,
2*(|alpha|^2+|beta|^2)=1 and |alpha|+|beta|<=1. Multilinearity
therefore preserves the largest basis budget for arbitrary normalized
TT leaves. Restoring the pair coupling gives

||H2_temporal||F<=530*(a+b)^2/(sqrt(kappa)*a*b).

Spatial rotation preserves the norm; stereographic endpoints are
handled by bounded limits. The estimate does not erase its explicit
energy-hierarchy factor, and gives no bound on a separate hard remainder.

## Three-ray planar family

Set the three half-angle coordinates to0,r,t, all in the xz plane,
and fix energies1/64,2/64,3/64. Compute the full four-tree
pure-three-graviton root and the complete degree-three nonlinear
pullback from notes/chart.md. All eight plus/cross choices are included.

For the total Q, exact algebra gives

Q^2=D/[1024*(1+r^2)*(1+t^2)],
D=8r^2-12rt+9t^2+5r^2t^2
 =2*(r^2+t^2)+6*(r-t)^2+t^2+5r^2t^2.

Thus D>=2*(r^2+t^2). After multiplying EACH temporal tensor entry by D,
its denominator factors only as a positive rational constant times
powers of1+r^2 and1+t^2. Every numerator monomial has total degree
at least two. These are exact symbolic factorization and support tests,
not fits, floating samples or assumptions.

For |r|,|t|<=1, every monomial |r|^i |t|^j with i+j>=2 is at most
r^2+t^2. If i>=2 or j>=2 this is immediate; the only remaining case
has i=j=1 and follows from2|rt|<=r^2+t^2. Therefore the numerator is
bounded by its coefficient l1 sum times(r^2+t^2).
The positive stereographic denominator factors are at least one, so
division by D gives an entry bound equal to coefficient_l1/(2*constant).

Summing all16 entry bounds (including both offdiagonals) bounds
the entrywise l1 norm and hence the Frobenius norm. For the ordered
basis triples000,001,010,011,100,101,110,111 the exact budgets are

25982722/225,122552/3,181128/5,828556/25,
5200136/225,4130957/180,2862328/135,1179848/225.

The largest is25982722/225<120000. The same unit-TT l1 argument
as above extends this maximum to all normalized tensor leaves by
trilinearity. This proves a UNIFORM bound on the specified planar box,
including simultaneous approaches to the origin and pair-collinear
loci. It is stronger than checking a few directional limits.

The original point API remains undefined on its exact internal-pole
loci. The canceled algebraic representative can have bounded,
direction-dependent all-parallel limits; boundedness does not assign
a unique value to an original0/0 expression.

## Energy and coupling homogeneity

A pure-gravity current with N free leaves and a propagated root has
V two-derivative Einstein vertices and exactly V internal propagators,
including the root propagation. Uniform momentum scaling contributes
lambda^(2V)*lambda^(-2V)=1.

If r_v is the valence at vertex v, canonical vertex factors multiply
to kappa^(sum_v(1-r_v/2)). A tree with N+1 external legs and V-1
unpropagated graph-internal edges has sum_v r_v=N+1+2(V-1).
The exponent is therefore-(N-1)/2. Propagators in canonical
normalization contribute no extra kappa. The projection Pi_Q is
homogeneous of degree zero in momenta. Hence N3 scales as1/kappa.

A generic exact common-energy/kappa rescaling test independently checks
these identities. The bound is independent of a common positive energy
scale but not of arbitrary energy ratios: the proof fixed1:2:3.

## Remaining problem

The positive polynomial decomposition above uses coplanarity and those
fixed energies. It is not an arbitrary three-dimensional estimate and
does not prove a bound uniform over soft energy hierarchies. Higher
multiplicity and complete hard-amplitude norms require further work.
No phase-space integration, virtual completion or inclusive probability
sum is supplied by either current bound.
