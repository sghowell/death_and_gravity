# S5.8.CD — Finite-time nonexceptional M1 tree control

Operational contract selected 2026-09-05, after the earlier reduction and
normalization results. It is not a preregistered observable and does not
change the frozen covariant witness or linear classification.

## Object and admitted modes

Use the fixed CD/M1 witness in its physical matter metric, the cubic and
quartic phase Hamiltonians of S5.6.CD, and the exact coupled canonical map
and free propagation theorem of S5.7.CD. At any finite t0 use fixed units
ell0=tau*sqrt(1+u0^2), normalize a(t0)=1, and take

    I0=[t0-ell0/100,t0+ell0/100],
    10^11 <= |k| <= 10^12.

Here k/ell0 is the initial physical momentum. Signed external spatial
momenta sum to zero and every nonempty proper subset sum has norm at
least 10^11. This is a **hard-transfer restriction on the observable**;
soft modes are not removed from the underlying action. A quartic internal
wavevector has initial norm at most 2*10^12. All angular configurations
meeting these restrictions are included, not merely example momenta.

Choose one chart throughout I0: gamma for |x0|<=9/50, unitary otherwise.
S5.7 proves its admissibility. The two normalized scalar coordinates obey

    Ydot=P-Omega*Y, Pdot=-W*Y-Omega*P,
    E=(|P|^2+Y^dagger*W*Y)/2.

At the left endpoint put R=W(left)^(1/2), the positive matrix square root,
and choose the two free mode columns

    U0=(2R)^(-1/2), P0=-i*(R/2)^(1/2).

Evolve with the **exact coupled free equations**. In particular
Ydot0=P0-Omega(left)*U0, not P0. Each tensor polarization uses the
corresponding scalar positive-frequency initial data. This specifies a
local Gaussian complex structure on the tested modes; it is not generally
the ground state of Hnorm=E-P^T*Omega*Y. No global vacuum, Hadamard
extension to all momenta, in/out state prescription, or WKB approximation
is claimed. Both scalar columns can mix both original scalar species.

## Observable and criterion

Work on R^3 with Fourier measure d^3k/(2*pi)^3 and the spatial gauge of
S5.6, with gauge transformations decaying at infinity. Decompose one- and
two-particle Hilbert spaces into fixed total spatial momentum fibers.
The criterion is operator norm at most 1/1000 for:

1. the connected cubic one-to-two and two-to-one transition blocks;
2. the connected quartic-order hard two-to-two block, retaining one H4
   insertion and two time-ordered H3 insertions.

Include both scalar and both tensor internal free-mode columns. Retain the
entire canonical inverse map, gamma swap, both mixed matter shifts, all
matter-sourced spatial constraints through degree three, and every
Hamiltonian boundary term that can reach cubic or quartic order.

Only the specified tree contributions are bounded. Disconnected spectators,
vacuum production, contractions within a vertex, loops, ordering-dependent
loop effects, and higher field orders are excluded. The hard-transfer mask
is part of block 2; no bound on its removed forward/soft entries is inferred.
No total-momentum delta distribution is estimated in absolute value.

The target is **one finite sufficient M*tau**, independent of center time,
meeting the criterion. An exact positive-series recurrence computes the
kernel and Schur bounds and the first decimal power meeting the inequalities,
M*tau=10^324.
It need not be a realistic or near-optimal scale; no background duration or
observational parameter value is fixed by this mathematical test.

## Acceptance and limits

The certificate must replay the coupled free-control input, verify the
canonical initial CCR, give explicit rational canonical-seed and invariant
coefficient bounds, include all matter/York terms in its finite recurrence,
and certify the final scale inequalities. Written coefficient induction,
free-ODE and continuum Schur proofs support the exact arithmetic.

This closes the stated M1 finite-time cubic/quartic tree gate only. It does
not close all-orders perturbative control, infrared/zero-mode problems,
inclusive scattering, a physical cutoff or optimized hierarchy, nonlinear
PDE/BKL/backreaction stability, radiative protection, UV matching or full
P8. No gravitating periodic box with unresolved homogeneous constraints is
used as a substitute for the continuum fiber argument.
