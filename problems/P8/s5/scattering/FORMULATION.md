# S5.3b operational contract, v1 — 2026-09-04

This fills S5.3's required frequency, IR/angular, time and order scope for
the fixed D-only M0 family. It does not change the witness functions or the
completed linear classification, and does not close unrestricted control.
The operational criterion was chosen during this continuation, not before
observing the earlier algebraic results.

At any cosmic time t0 use fixed units `ell0=tau sqrt(1+u0^2)` and normalize
`a(t0)=1`. Work on R^3 with spatial gauge transformations decaying at
infinity, and Fourier measure `d^3k/(2*pi)^3`; k/ell0 is the initial physical
momentum. The window is `I=[t0-ell0/100,t0+ell0/100]`. Decompose the one-
and two-particle Hilbert spaces as direct integrals at fixed total spatial
momentum P. No physical periodic universe or zero-mode charge constraint
is discarded to define these states.

Allowed external modes have `10^8<=|k|<=10^9`. Signed spatial momenta sum
to zero and every nonempty proper subset sum has norm at least `10^8`.
The single internal wavevector of a four-leg tree has norm at most `2*10^9`.
This hard-transfer restriction concerns the observable, not removal of soft
modes from the action. Forward/zero-mode and loop effects are not bounded.

Choose gamma canonical variables for `|x0|<=9/50`, unitary otherwise, and
the exact volume-normalized oscillator variables Y,P of S5.3a. The chart
remains admissible over I. At the left endpoint choose free initial data
`u=1/sqrt(2 Omega)`, `udot=-i Omega u`, then evolve with the **exact** free
oscillator equations. This defines local free Fock blocks, not a global
vacuum, in/out S matrix or chart-independent frozen-state prescription.

The criterion requires norm <=`1/1000` for (i) connected cubic one-to-two
and two-to-one particle blocks, and (ii) the connected quartic-order hard
two-particle transition block, uniformly in the allowed total-momentum fibers.
For (ii) retain one H4 insertion and two H3 insertions, with scalar and both
tensor internal polarizations. Exclude disconnected spectators and loops.
The hard-transfer mask is part of block (ii); no bound on removed entries
is inferred. Vacuum production, the full Fock-space norm and quantum-ordering
differences first contributing loops are outside this tree criterion.

Bound the ordinary kernels in each momentum fiber by a finite-measure Schur
estimate. The total-momentum delta function is handled by direct-integral
decomposition, never by an absolute bound on a distribution.
The target is one common finite M tau meeting this criterion at all t0,
which the verified scale restoration permits. The sufficient value is not
an optimized physical estimate. Higher tree orders, loops and UV
admissibility remain separate gates.

Revision before certification: replace the exploratory physical-box/state-count
wording by momentum-fiber Schur bounds. On a closed gravitating torus, residual
homogeneous constraints require separate treatment; a box with those constraints
silently dropped would not establish a physical hard-channel theorem. The
revised criterion also excludes infinite-volume vacuum production explicitly.
