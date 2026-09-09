# Global signed wavefront condition

Write k=|xi| and use only k above the fixed initial cutoff.
The remaining compact-frequency sector is smooth after the exact
global evolution. All statements below concern classical symbols
with every u and xi derivative bounded on each specified compact
time strip. Spatial translation invariance makes symbol composition
ordinary matrix multiplication; no omitted spatial symbol terms occur.

## Local signed phase representation

On the initial canonical slab the preceding Riccati proof constructs
an exact negative-graph solution modulo a smoothing symbol. Its
configuration generator has principal part -i*k*V omega V^-1.
After the smooth V change it is -i*k*omega plus an order-zero symbol.
The two frequencies differ on every compact slab.

A formal change I+sum k^-j T_j diagonalizes this two-mode equation
modulo S^-infinity. At each order the known off-diagonal remainder
is divided by -i*(omega_i-omega_j); the diagonal remainder is retained
as amplitude transport, not set to zero. The explicit gap bounds in
domains.py justify each division with all time derivatives.
The same cutoff construction used for the Riccati series sums this
formal change into an invertible order-zero symbol for large k.

Each scalar diagonal evolution consequently has form

    exp(-i*k*integral_0^t omega_j(u) du) * a_j(t,k),

where a_j is an order-zero symbol with a full asymptotic expansion.
Indeed the diagonal remainder is an order-zero symbol; integrating
and exponentiating it leaves a classical order-zero amplitude on
a compact interval. Smoothing errors in the transformed equation
remain smoothing after Duhamel and the complete polynomial bounds.
Multiplying by the graph, initial covariance Gram factor, and original
density maps only changes finite amplitude orders. Therefore the
exact negative modes have the stated signed oscillatory forms
modulo smooth kernels on this slab.

## Across every actual chart transition

S6.98 gives the complete scaled velocity systems Y=(k*x,x') in
both charts with generators k*J0+R, where R is a classical order-zero
symbol at large REAL k. Its rational denominators are separated
uniformly; neither the gamma velocity pole nor a chart zero is crossed
using an invalid inverse.

The new exact overlap calculations are essential. The map from
canonical packets to gamma velocity variables and its inverse have
order zero, with leading map a^-3/2*diag(I,K^-1).
The gamma-to-unitary velocity map and inverse also have order zero,
with leading map

    diag(-Theta/Lambda,1,-Theta/Lambda,1).

At u=+-3/16 both pivots are nonzero. These leading maps preserve
each sign of the same physical clock and matter frequencies.
The full rational maps, not just these limits, are checked against
their exact inverses. They intertwine the complete evolution because
both are maps from the same original density generator.

For completeness, sign preservation to ALL orders follows from a
formal-projector uniqueness argument. In a balanced chart seek
P~P0+sum k^-j P_j with P^2=P and P'-[M,P]=0 modulo smoothing,
where P0 selects the two negative-time-frequency modes.
The formal diagonalization supplies such a projector. On the initial
slab the Riccati graph and its conjugate are complementary, because
R0-conjugate(R0)=-2i B0 is invertible; their algebraic projector supplies
the same formal solution.

If two such projectors first differ at order k^-n, the leading
invariance equation says that difference commutes with J0.
The opposite-sign gaps make its off-sign blocks zero. The
idempotence equation P0 D+D P0=D then makes its two diagonal-sign
blocks zero. Thus the difference is zero, a contradiction.
Induction gives uniqueness modulo S^-infinity. This argument does
not discard a nonzero order-zero transport term.

Conjugation by an exact order-zero overlap preserves both equations
and the leading sign projector. Uniqueness therefore proves that
the two chart constructions agree modulo smoothing. Differentiating
U(s,t) P(t) U(t,s) shows that the exact evolution transports this
projector modulo smoothing; the defect integral stays smoothing by
the full polynomial transfer bounds. Starting with the initial
negative graph, no nonsmoothing opposite-sign component appears
on the next segment. Apply this at the fixed overlaps. Every finite
interval uses at most three regular segments. Thus negative-frequency
microlocal solutions extend throughout the global clock, in both
time directions, without requiring an all-time uniform gap.

Within each sign block the two distinct frequencies can again be
diagonalized as above. The universal four-mode first-order identity
is checked in recursion.py; further orders follow the same
off-diagonal division and retained diagonal transport.
All such divisions are finite on every chosen compact strip by
the positive-polynomial gap bounds.

## Two-point wavefront set

On each finite pair of time neighborhoods the nonsmoothing part of
W is a finite sum of oscillatory integrals with finite-order classical
amplitudes and phases locally of the form

    (x-y).xi - |xi|*integral^t omega_j
                     + |xi|*integral^s omega_l
                     + a homogeneous spatial-momentum phase
                       from intermediate signed segments.

The exact value of the last phase is irrelevant to the sign
inclusion. First-time differentiation gives tau=-omega_j(t)|xi|;
second-time differentiation gives sigma=+omega_l(s)|xi|.
The spatial covectors are xi and -xi. Integration by parts away
from the stationary covectors gives rapid decay in every closed
cone disjoint from these sets, to arbitrary order because the
amplitudes have all symbol derivatives. This proves directly

    WF(W) subset V+ x V-,

with V+ as in the formulation. No component with only one zero
covector is left: nonzero xi occurs in both legs of the same
oscillatory integral; the compact-frequency part is smooth.
Finite-order original density maps do not enlarge this sign
inclusion.

The sets V+ and V- are conic, relatively closed in the punctured
cotangent bundle and disjoint. On compact strips all omega_j are
strictly positive. A sequence with xi tending to zero also has
tau tending to zero, so has no nonzero missing boundary covector.
The two characteristic sheets may approach each other at infinite
time without violating this local closedness.

Finally W-W^transpose=i*hbar E_Q/kappa. Its transpose has the
opposite wavefront signs, so this equality implies

    WF(E_Q) subset (V+ x V-) union (V- x V+).

Together with notes/green.md and positivity, this establishes the
generalized Hadamard assertion. We do not claim an exact
bicharacteristic relation or a one-metric Klein--Gordon parametrix.

Terminology follows [Fewster, Definition 5.2](https://eprints.whiterose.ac.uk/id/eprint/234321/1/Hadamard_arXiv_V3.pdf):
decomposability and a separated-sign two-point wavefront condition
define the generalized state class. The normally-hyperbolic
existence theorem in that paper is NOT used; this construction
supplies the state and the needed decomposability itself.
