# P8 S6.78: seven-mode finite-time hard-tree control

Research checkpoint, 2026-09-08. Original P8 is OPEN.
This uses the same classical seven-mode action and fixed
margin as S6.77. It adds no action deformation, state
reselection, order reduction or higher-operator prescription.

## Observable and fixed domain

Use any interval of full length 1/(2*10^9) contained in
I=[-1/2,1/2]. Its half-width is 1/(4*10^9). In constant
spatial coordinates with a(center)=1, external momenta
have magnitudes in [10^12,10^13]. Every nonempty proper
transfer occurring in the cubic/quartic tree is at least
10^12. The internal exchange momentum is at most 2*10^13.

The observable comprises the finite-time cubic 1-to-2 and
2-to-1 blocks and the hard-masked connected 2-to-2 tree
block through H4 and H3 squared. Use d^3k/(2*pi)^3 and
fixed-total-momentum fibers on R^3. All seven physical
mode columns, including massive Proca modes, are retained.

Scalar columns have canonical initial data defined from
a real Cholesky root of their actual positive kinetic
matrix and its symmetric momentum boundary. Their
antisymmetric initial-velocity contribution is retained.
Tensor columns use the declared canonical oscillator
data. The vector columns remain the original selected
Borel-prepared Gaussian modes; their nonzero mixing
with the WKB comparison is not set to zero.

The scalar/tensor initial data specify finite-band
canonical columns, not a new all-momentum Hadamard
or semiclassical state theorem. The block construction
uses mode-kernel Wick contractions. It does not require
or assert globally implementable homogeneous free
evolution on an infinite-volume Fock space.

## Quantitative statement

The full two-scalar time-dependent Lagrangian keeps
the antisymmetric mixing matrix. Its principal
gradient and kinetic lower bounds, finite-q remainder
and exact energy identity give a uniform free-column
bound on both charts. The tensor and selected-state
vector bounds complete the seven-mode phase estimate.

A common component seed 10^42 bounds the normalized
phase columns in the raw variables of S6.77, before
their overall 1/(M*tau) factor. Keep both scalar-matter
boundary shifts, the gamma canonical swap, the
longitudinal Proca normalization, tensor duals and the
constant-center Fourier change. The fixed-coordinate
Hamiltonian also retains (a/a_center)^3<=8.

The positive H3/H4 recursion, Wick/species counts,
both time orderings and a fixed-fiber Schur bound
give explicit rational constants C3 and C4:

cubic block norm <= C3/(M*tau),
connected quartic tree block norm <= C4/(M*tau)^2.

The separately named scale example S6.78-P353,
M*tau=10^353, makes each bound at most 1/1000.
This is a conservative sufficient example, not an
optimal or necessary hierarchy.

## Boundaries

The band lies above m*tau=1000 and retains the heavy
vector modes. It is not a domain in which the vector
has been integrated out. No light-only heavy-mass
expansion is used there.

The old M*tau=10^24 quantum-response and frozen-pole
examples are not replaced or transferred to P353.
The selected dimensionless vector state and mass
remain fixed; only this named sufficient scale
example is different.

No full Fock-space norm, disconnected vacuum/spectator
contribution, loop graph, zero/forward channel,
infinite-time scattering amplitude, uniform nonlinear
spatial existence, Wilsonian/all-orders cutoff or
UV completion is claimed.

Finite Wilson matching, omitted-operator/loop/threshold
budgets, common-parent V/G/B and original P8 closure
remain open.
