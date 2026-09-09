# Cauchy propagation, the Proca correction and both signed cones

The physical and auxiliary spacetimes have smooth positive
scales a,b>=1. Any causal curve parameterized by u obeys
|dx/du|<=1. A finite terminal u therefore forces a finite
spatial limit and the curve can be extended in the smooth
metric. Thus an inextendible causal curve traverses the
entire real u-axis, and each R^3 constant-u surface is Cauchy.
The same argument applies to a sufficiently small common
slab about u=0.

For the Proca sector, apply the Cauchy-slab extension and
Hadamard propagation statement in
[Moretti, Murro and Volpe, Proposition 4.7](https://arxiv.org/pdf/2210.09278),
with the following explicit repair. Its use of the Proca
projection as an elliptic operator is invalid: the symbol is
nonzero but has a kernel. The needed Proca causal-propagator
wavefront equality is established by
[Fewster, Theorem 5.1](https://arxiv.org/pdf/2503.12544),
using polarization information. The repaired theorem, not the
unrepaired ellipticity argument, is used here.

The native negative control takes a null covector n=(1,0,0,1).
For eta=diag(-1,1,1,1), the symbol n(n^T eta)/m^2 has rank one,
zero determinant and square zero. It is emphatically not
invertible. Merely checking a nonzero entry would not justify
an elliptic regularity inference.

The neutral Proca operator has the same constant positive
mass and positive-action normalization as our Hamiltonian.
The reference convention g_F=-g gives delta_F=-delta here:
its -delta_F d+m^2 is our delta d+m^2. This flips metric
conventions, not the physical theory. Time orientation is
fixed by the positive-energy flat covariance with e^-i omega t,
then propagated continuously; no future/past label is imported
without this convention check.

The flat covariance's explicit three-polarization factor is
the standard positive Minkowski massive Proca two-point
distribution in this normalization. It is Hadamard on the
flat reference past. First propagate across b, then restrict
to the shared slab and propagate with a. Proposition 4.7(a)
also fixes uniqueness on each Cauchy extension. This proves
the state prescription of the previous note has the
physical metric's signed Hadamard relation globally.

For the minimal scalar used to realize the tensors, use
[Fewster, Theorem 5.4(c)](https://arxiv.org/pdf/2503.12537)
for Hadamard propagation from the flat reference past
through the same two Cauchy extensions.

Microlocally near the nonzero massless characteristic cone,
the spatial TT projector is an order-zero smooth homogeneous
symbol. Since nonzero characteristic covectors have
nonzero spatial momentum there, the spatial operator may
be treated as an order-zero spacetime pseudodifferential
operator on that conic neighborhood. It cannot introduce
another wavefront cone. The infrared piece is smooth by
the independent integral argument, and the rank-two
principal projector does not remove the entire matter
null relation. Thus the physical tensor distribution
retains the matter cone.

The scalar factor's two-cone Hadamard result is inherited
without changing its Borel or infrared choice. The added
sector cross-covariances vanish. For a finite block-valued
distribution, the wavefront set is the union over its
blocks: smoothness of the matrix is equivalent to smoothness
of every component. Hence the full phase has exactly the
union of the inherited clock and physical matter signed
relations. Positive and negative first-leg frequency
directions remain disjoint because u is a common global
time. A constant massive pole adds no principal cone.

This is a direct blockwise proof for the declared physical
phase. It does not assert that the full constrained gravity
system belongs to a formal Green-hyperbolic operator category
without further work, nor that nonlocal TT gauge components
are local observables at spacelike separation.
