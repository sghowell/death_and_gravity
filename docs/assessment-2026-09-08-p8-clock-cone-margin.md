# P8: background-preserving clock-cone margin

Original P8 remains open. This follows the
[controlled vector-state audit](assessment-2026-09-08-p8-hadamard-vector-state.md).
No frozen action certificate or original closure contract changes.

## New named action

[S6.56](../problems/P8/s6/matching/affine/kinetic/aligned/margin/FORMULATION.md)
adds the recorded physical lower scalar
Delta P=(M²/tau²)*epsilon*(x+1)²/h², with
h=(1+u²)³ and 0<epsilon<=1/100.

Its value and first variations vanish on x=-1, preserving the
exact classical bounce and free-matter history. The full
volume variation changes the lapse-square coefficient by
delta_J=4epsilon/h². The new regular Hamiltonian and full
rolling auxiliary Hessian are recomputed. The local nonlinear
degree count persists, without a nonlinear stability claim.

The two scalar squared principal speeds are exactly
1 and J/(J+delta_J). Both the punctured unitary and regular
gamma charts are derived with their moving boundary terms.
They cover the bounce; the clock is positive and strictly
subluminal at every finite time. Matter remains luminal.

Positive polynomial coefficients prove J*h²>=1199/800
globally and J*h²<36 on |u|<=1/2. For epsilon=10^-6,
the fractional added light kinetic form is below three parts
per million and the compact clock-cone margin exceeds 10^-7.

The fixed finite vector potential has B_xx=2392/(6561h²),
hence negative Lorentzian P_xx. At the same L=10^12,R=1000
example it is below one billionth of the added positive
P_xx in magnitude. This is an isolated local jet comparison,
not a full constrained quantum cone. No finite tadpoles or
counterterms have been cancelled.

The vector clock operator, state, matching and first-variation
bounds are unchanged. Old light interaction/response results
are not automatically transferred to this changed action.

## Remaining work

The still-luminal directions have no strict margin against
arbitrary unsigned corrections. Quantum feedback, corrected
background/constraints/cones, interacting cutoff and V/G/B
remain open. A first-order homogeneous response to the
fixed bounded vector stress is now derived and undergoing
its separate certificate checks; it is not a self-consistent
quantum solution. No user-intervention blocker is present.

## Verification

The report pins 13 sources and fully rebuilds S6.55.
It checks 33 named identities comprising 62 scalar entries,
15 continuous proof checks and 51 rejected inputs.
Focused science passes 9 tests in 2.51 seconds; the ordinary
suite passes 43 tests in 440.39 seconds without the broad
GCD adapter. The separate read-only CLI passes.

The full P8 regression through S6.56 passes 5008 tests in
956.19 seconds, without excluding a frozen checkpoint. The
exact GCD adapter passes 128 original tuple comparisons and
records 5879 domain fallbacks and 6509 exact descents. This
run collected before the later S6.57 working tests were added.

The report SHA-256 is
b6952b45787005033a7ae19176c4668525044b37a2afb8753651ca99915d47dd.
The evidence consists of exact algebra and written continuous
proofs, not proof-assistant formalization or original P8 closure.
