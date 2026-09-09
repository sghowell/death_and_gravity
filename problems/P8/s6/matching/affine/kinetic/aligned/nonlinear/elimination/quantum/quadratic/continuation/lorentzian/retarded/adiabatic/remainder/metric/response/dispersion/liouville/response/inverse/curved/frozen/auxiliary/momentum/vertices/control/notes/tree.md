# Fixed-fiber, finite-time hard connected trees

## Normalization and the retained Hamiltonian

Use the complete classical reduced H3 and H4 of S6.77.
The chosen canonical columns solve the full actual free
two-scalar, two-tensor and three-vector systems. A
linear time-dependent canonical change adds only a
quadratic generator; it does not delete or invent a
higher interaction vertex.

Set L=M*tau. In dimensionless units the action has
overall factor L^2. Rescale both members of every raw
phase pair by 1/L. For a homogeneous phase polynomial,

L^2*H_n(Q/L,p/L)=L^(2-n)*H_n(Q,p).

This full-phase homothety, including momentum-dependent
vertices, is checked exactly at n=2,3,4. The physical
Hamiltonian's 1/tau cancels dt=tau*du in its time
integral. No dimensional cutoff is inferred from it.

The raw majorant is evaluated at S=10^42, K=8*10^13,
D=2*10^-12, with all vector and matter terms enabled.
It bounds every labelled complex Fourier H3/H4 kernel
in the stated transfer domain, including factorial
multiplicities. Let these bounds be h3,h4.

The Hamiltonian is integrated in constant-center
coordinates. Although each kernel is expressed in local
physical components, its density must be multiplied by

(a/a_center)^3<=8.

This factor is distinct from the square-root Fourier
normalization of each individual phase column. Both
are retained. Define B3=8*720*h3 and B4=8*720*h4.

## Combinatorics and scope of the tree kernel

The factor 720=6! conservatively covers vertex label
permutations, few-particle symmetrization and Wick
assignments. A cubic or quartic labelled coefficient
already contains its appropriate factorial, so this
is an overcount, not compensation for a missing 1/n!.
It is valid also when physical species labels repeat.
For example, a cubic step from two to at most five
quanta has occupation factor at most sqrt(5!/2!) and
at most 3! assignments: 6*sqrt(60)<720. A quartic
two-to-two step has bound 4!*2=48<720. Keeping all
species equal maximizes these occupation factors.

For a connected four-external-leg two-cubic tree, choose
one of three unordered pair partitions, one of seven
internal physical mode columns, either internal complex
orientation and either temporal placement of the two
vertices. Use the separate bound 3*7*2*2=84. Each endpoint
of the internal contraction is one of the normalized
phase columns bounded by S. The internal momentum is
fixed by the external pair; there is no independent
internal momentum integral for this connected tree.
The full square of side T bounds either time-ordered
triangle and the product of its two kernel bounds.

This mode-kernel proof does not restrict intermediate
particle number to one or three: counterrotating cubic
terms can pass through five-particle states. Connected
Wick contraction accounts for them without pretending
that the truncated particle subspace is invariant.
Vacuum loops and disconnected spectators are outside
the observable, not silently included in its norm.

Contract canonical Hamiltonian phase fields directly.
In particular, p-p contractions are ordinary products
of exact free phase columns. Do not replace them with
time derivatives of a Feynman coordinate propagator
and then omit its distributional contact term. The
full canonical H4 already includes the actual lapse
response and stationary/Legendre corrections; adding
an ad hoc Lagrangian contact a second time is incorrect.

## A deliberately generous fixed-fiber Schur bound

Use d^3k/(2*pi)^3 and decompose by conserved total
spatial momentum P. The one-particle fiber is C^7.
The two-particle fiber has one relative momentum and
two species labels, with bosonic symmetry if desired.
Restriction to the symmetric subspace cannot increase
the norm. All supported relative momenta fit the box
[-2U0,2U0]^3, including exchange momenta. Dropping
(2*pi)^-3<1 only increases the measure bound.

Write mu=7*(4U0+1)^3>7. The finite species/momentum
measure of each two-particle fiber is at most 7*mu;
that of the one-particle fiber is seven. A kernel
bounded uniformly by B therefore has norm at most
B*sqrt(measure_input*measure_output), by Cauchy-Schwarz,
or at most B*max(measures), by Schur's test. Both are
bounded by B*Schur, with the convenient common

Schur=8*mu^3.

Hard external and transfer masks have magnitude at
most one and cannot increase these absolute-kernel
estimates. The bounds are uniform in P, including empty
fibers. There is no infinite-volume delta(0) norm claim.

For T=1/(2*10^9), obtain

C3=Schur*T*B3,
C4=Schur*(T*B4+84*T^2*B3^2),

||cubic 1-to-2 or 2-to-1 block||<=C3/L,
||hard connected 2-to-2 tree block||<=C4/L^2.

The stored constants are exact rationals, not floating
estimates. Increasing L through powers of ten gives
the sufficient named example L=10^353, at which both
bounds are at most 1/1000. This search minimizes only
the decimal power satisfying these particular loose
majorants. It does not find a physical strong-coupling
scale, an optimal hierarchy or a necessary condition.

## What this does not close

The construction uses finite-band canonical field
columns and Wick kernels, not a claim that homogeneous
evolution is globally unitarily implementable on an
infinite-volume Fock space. It is not an infinite-time
S matrix or a full Fock-space operator-norm bound.

The old L=10^24 response and frozen-pole examples remain
their own examples. The new P353 hierarchy is explicitly
named, keeps the same dimensionless action, margin,
vector mass and selected vector state, and supplies
only the stated classical finite-order tree result.

Even a very small raw field amplitude does not prove
uniform nonlinear York solvability at higher orders;
new zero subset transfers and spatial analytic-domain
issues must be addressed rather than extrapolated away.
No higher-loop, higher-operator, threshold matching or
all-orders error budget is established here. In
particular this is not a Wilsonian cutoff or a
quantum-stability/cone claim. Common-parent V/G/B and
the original P8 problem remain OPEN.
