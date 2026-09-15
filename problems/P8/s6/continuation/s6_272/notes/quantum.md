# Full operator bounds and a finite volume turnaround

## The new volume amplitude, with all ordering contacts

The full complex volume amplitude is A=10^-320 on the initial ball 4R.
Keep the same original covariance whitening, two explicit cutoffs and
96 phase coordinates. The S6.271 high-phase proof applies without
changing its geometry or derivative constants:

    J_n=2A 2056^n(n!)^3/R^n, 0<=n<=196.

All 196 successive ratios are less than 10^-9. This is a phase
derivative statement, not an assumption of 196 time derivatives of
the original C5 profiles. The simultaneous Cauchy radius is R/8,
not R/4 in all 96 coordinates.

For the full calibrated-coherent volume,
||F_C-I||<=A+(96/4)J_2<2A. The complete heat remainder gives

    ||F_W-F_C|| <=10^112 (96^2/32)J_4.

The explicit normalized Gaussian-frame Schur constant, all mixed
derivatives, cutoff contacts and implicit-source derivatives are
included. Thus ||F_W-I||<epsilon=10^-264 and F_W>I/2.
This is a concrete operator estimate, not symbol supremum being
identified with operator norm or generic positive Weyl quantization.

## Exact unitary evolution without small-displacement claims

For either ordering, the entire real extended Hamiltonian is a
retained scalar plus a bounded Schwartz-kernel operator, with
norm-continuous time dependence. The full bounded Dyson construction
therefore gives the exact unitary, including its scalar phase.
The inherited seminorm argument gives the common Schwartz strong
equation. The complete source smoothness also makes these readout
means C1. Only qualitative finite time-derivative bounds are needed
for that regularity, not a small acceleration error estimate.

The available integrated Hamiltonian norm bounds on this interval
are much larger than two. Hence the old tiny-time unitary, evolved
leakage, ordering-state and cutoff-state estimates do NOT extend here.
The proof below uses only exact unitarity and a uniform operator
bound. Independent noncommuting state-flip fixtures explicitly show
why a state need not stay near its initial vector.

Use the original metaplectic reference consistently: physical states
are U_ref U_I psi0 and the physical readout is U_ref F U_ref^-1.
The expectation equals the interaction-picture one. The unchanged
full covariance, scalar phases and reducing three-translation-charge
sector remain; no coherent-label conditioning is performed.

## Endpoint comparison and interior minima

For each ordering and cutoff, every corresponding normalized evolved
state satisfies

    1-epsilon <= <F(u)> <=1+epsilon.

The unchanged physical reference volume is a^3=(1+u^2)^6. Therefore

    nu(0) <=1+epsilon,
    nu(+-T) >=(1+T^2)^6(1-epsilon),
    nu(+-T)-nu(0)>5T^2=5 x 10^-260.

A continuous mean on the compact interval has a minimum, and neither
endpoint can be a minimizer. If u_star is any global minimizer,
nu(u_star)<=nu(0), so

    6(1-epsilon) u_star^2 <=2epsilon.

Thus every minimizer lies strictly inside |u_star|<10^-132.
The C1 mean is critical there, and the mean-value theorem supplies
a point of negative derivative before it and a point of positive
derivative after it. It does not prove monotonicity on either whole
side, uniqueness, or strict positive acceleration at a minimum.

This is a theorem about the declared finite, phase-extended readout
with EXTERNAL homogeneous functions. It does not show that the state
remains in the classical core on this longer interval, nor that the
homogeneous functions solve quantum mean equations. Original
unlocalized dynamics, regulator removal, physical matching,
omitted loops, UV/Regge control and P8 closure remain open.
