# Exact bounded evolution, evolved leakage and two-regulator readouts

## Same seed and entire propagator

For each explicit cutoff let A_c(u)=Q_V[(1-D)g_ext,c(u)].
The common full phase-domain proof and cutoff estimates imply a bounded
self-adjoint, norm-continuous generator with ||A_c||<B=10^1010 on
real |u|<=T=10^-2000. Its actual scalar background phase is retained.

The S268 norm-convergent Dyson construction therefore yields the
exact unitary interaction propagator U_I,c. Each term has norm at
most (B|u|)^n/n!, and the full integral equation with unitarity gives
||U_I,c(u)-I||<=B|u|<=10^-990.
This is a bound on the entire unitary, not an evolution defined by a
first-order or finite occupation truncation. The complete selected free
metaplectic propagator gives U_full,c=U_ref U_I,c. The inherited
Schwartz-domain argument applies to the smooth compact coherent part
plus the retained scalar identity, and supplies the strong equation
on that domain for the declared regulated Hamiltonian.

The full symbols, same-seed coherent transform and whitened radial
cutoffs commute with all three original translations. Their common
zero-charge factor is reducing under the exact propagator. The seed is
the same actual translation-invariant pure preparation, not a
conditioned or projected Gaussian. Coherent integration labels need
not themselves have zero classical charge; that is not an omission of
the original residual quantum constraint.

## The initial tail is positive, and its evolved bound is not zero

In the SAME pure-state whitened coordinates the Husimi probability
density is the standard96-dimensional normal. The core has R=10^20
and d=48, so t=R^2/2 and its exact tail is

tail=exp(-t) sum_(k=0)^47 t^k/k! >0.

The inherited exact Chernoff calculation gives tail<exp(-10^39).
Because e^3>10 (already the first four positive exponential terms
sum to13), log10<3 and10^39>12000, this implies tail<10^-4000.
This rational surrogate is conservative and is never rounded to zero.

Let E_out=Q_V(1_outside-core), a positive contraction. For normalized
vectors, an expectation of a contraction changes by at most twice
the vector norm difference. Therefore

< U_I,c psi0, E_out U_I,c psi0 >
<=tail+2||U_I,c psi0-psi0||
<10^-4000+2*10^-990<10^-980.

The estimate concerns the evolved POVM probability. It does not claim
joint sharp support for noncommuting observables, an invariant classical
ball, a state projection or a uniform dimension/volume limit.

## Both full generators, all cutoff derivatives

For c=1,2, delta a=(1-D)[(chi_1-chi_2)(g-g0)].
It vanishes on the common core, including every derivative contribution:
the two smooth cutoffs coincide with the same constant1 there and are
flat at the core boundary. Their common scalar background cancels
in the difference, while remaining in each full evolution.
The complete symbol bound is ||delta a||<=2B.

The coherent Schwarz inequality gives
||Q(delta a) psi0||<=2B sqrt(tail).
Let one evolution be the comparison propagator. Duhamel and its
unitarity, followed by ||U_I,2(s)psi0-psi0||<=B|s|, yield for both
time directions

||U_I,1(u)psi0-U_I,2(u)psi0||
<=2B sqrt(tail)|u|+B^2 u^2.

At T the displayed exact rational upper bound is<10^-1970.
This compares two explicitly defined regular Hamiltonians; no
undefined original singular Hamiltonian occurs on either side.

## Compare the readout as well as the state

Consider either the two uncalibrated positive coherent volume operators
or the two separately calibrated volume operators. All have norm<2.
Their difference symbol vanishes on the common core and has norm
<=4E, E=10^-255; this bound also includes the calibrated derivative
difference, since each calibrated symbol is within2E of1.

For the same initial seed the absolute mean difference of this symbol
is<=4E tail. In a state evolved by the second regulator it is at most
4E tail+8E B|u|. The remaining change of state in the first
operator's expectation is<=4 times the preceding state difference.
Combining both pieces gives

mean_difference<=4 state_difference+4E tail+8E B|u|<10^-1230.

This compares the same type of readout across the two regulators;
it does not equate the positive and calibrated readouts to each other.
The result includes evolved operator-symbol difference, not just the
initial tail or a state-distance estimate.

## Physical picture and limits

The physical-picture observable is U_ref Q(F_ext) U_ref^dagger,
or the corresponding separately calibrated operator. The physical
state is U_ref U_I psi0. Unitarity cancels both reference factors in
the expectation, giving exactly the interaction-picture formulas above.
No fixed-window metaplectic covariance is invoked without also
transporting that window.

The finite-dimensional noncommuting matrix fixtures in the tests check
Duhamel, phase retention, leakage, readout change and physical-picture
conjugation. They are independent diagnostics, not a substitute for
the infinite-dimensional CCR/coherent-map theorem and not a numerical
simulation of the full P8 evolution.

Neither a long-time result nor regulator removal follows. B is large
and T deliberately tiny. Finite-time control of two regularizations
does not supply original Wilsonian matching, all omitted loops,
backreaction or nonlinear global completeness.

