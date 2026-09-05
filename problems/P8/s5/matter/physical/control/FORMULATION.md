# S5.7 — Coupled CD/M1 normalization and local free-energy control

Frozen scoped target, 2026-09-05. This is a new, separately pinned extension
of S5.6.CD; no prior source, certificate, or claim is revised in place.

## Input and question

Use the fixed CD/M1 covariant witness and **physical matter metric** of
P8-2.CD, with the full momentum-reduced quadratic Hamiltonian of S5.6.CD.
Keep both scalar species, their finite-momentum mixing, the physical volume
measure, and every time-dependent canonical generator. The principal
identity K=G is not a formula for the full finite-momentum oscillator.

The question here is whether exact coupled canonical coordinates and a
uniform sufficient local free-propagation band exist at every background
time. It is not yet an interacting strong-coupling or UV question.

## Acceptance contract

For u=t/tau, a=(1+u^2)^2, ell=tau*sqrt(1+u^2), and fixed nonzero comoving
momentum, set x=u/sqrt(1+u^2), q=ell^2*k_com^2/a^2, z=1/q. The endpoints
x=+-1 are compact limits, not points at finite physical time.

1. Derive the exact two-scalar canonical map from the full Hamiltonian,
   including the volume generator, symmetric momentum boundary, and
   antisymmetric connection. No instantaneous eigenbasis, WKB prescription,
   fixed-local-q differentiation, or deletion of mixing is permitted.
2. Exhibit covering positive kinetic charts: unitary for abs(x)>=1/9 and
   gamma for abs(x)<=1/4, both with q>=1000. The gamma condition remains
   q*Lambda^2>J+w^2/2. A phase chart and its velocity chart are distinct.
3. Write the normalized free Lagrangian as

       L2 = |Ydot+Omega*Y|^2/2 - Y^dagger*W*Y/2,
       Omega^T=-Omega, W=(q*I+M)/ell^2.

   Establish exact, uniform coefficient majorants for M and its covariant
   derivative. The same contract includes each tensor polarization.
4. Give a sufficient q band with positive W and a quantitative free-energy
   estimate on a specified local physical-time interval. A finite sampling
   grid or the principal cone alone does not meet this requirement.
5. Independently recover finite-q bounce formulas and nonzero-time fixtures
   by differentiating the original cosmic-time Hamiltonian, before point
   evaluation. Retain explicit omission controls and a low-q control that
   separates positive kinetic energy from positive oscillator potential.

## Result and scope

The replay proves the conservative sufficient threshold q_star=10^20. At
any center t0, take q0>=2*q_star and the window abs(t-t0)<=ell0/100. Choose
one chart for the whole window: gamma if abs(x0)<=9/50, unitary otherwise.
Then

    q/(2*ell^2) I <= W <= 3*q/(2*ell^2) I,
    abs(E_dot) <= 18*E/ell,
    7/11 <= E(t)/E(s) <= 11/7 < 2,
    E=(|P|^2+Y^dagger*W*Y)/2, P=Ydot+Omega*Y,

for every nonzero exact free solution and all s,t in that window. The
energy is chart-specific; no equality of energies across chart changes is
claimed. The threshold is an intentionally loose coefficient-majorant
bound, not a fitted, optimal, or physically recommended momentum scale.

This settles the coupled quadratic-normalization/free-propagation gate.
It does **not** settle normalized M1 cubic/quartic majorants, finite-time
interacting tree control, a cutoff hierarchy at a chosen M*tau, scattering,
loops, all-wavelength or nonlinear stability, UV matching, or full P8.
No fixed comoving mode stays in the stated high-q band for all time: q
falls like (1+u^2)^-3 along either tail.

The certificate combines exact algebra and rational coefficient bounds
with the explicit analytic proof in `notes/normalization.md`. Its proof
level is CERTIFIED in the repository's scoped sense, not Lean FORMALIZED.
