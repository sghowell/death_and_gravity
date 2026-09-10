# Formal source identity and stable-light LSZ

At any finite perturbative order, change variables in the regulated
generating integral using Phi=F(Psi). The exponent contains the full
S[F(Psi),H], the original counterterms after substitution and the
physical source J F(Psi), not J Psi. Retain det(F') until the
regularized perturbative argument is made.

An infinitesimal version follows by integrating the total functional
derivative of R times the source-weighted integrand. Its terms are
the variation of the action, divergence of R (the Jacobian) and
variation of the observable or source. Iteration gives the finite
formal identity. It is an order-by-order identity, not a
nonperturbative measure construction for a derivative map.

Here R' is a finite local differential polynomial. Expanding
Tr log(1+R') at any finite order gives ghost loops with identity
propagator and polynomial loop momenta. Each loop is scaleless in
dimensional regularization and vanishes, even with fixed external
momenta as polynomial coefficients. There is no inverse heavy
kernel in this Jacobian: R is local. This is done before removing
the regulator. It does not license dropping generated matter
interactions, transformed sources or composite counterterms.

The native independent diagnostic uses the ordinary Gaussian
integral, where the Jacobian does NOT vanish. For F(y)=y+r y^3,
its measure includes F'=1+3r y^2 and

    exp[-(S2(F)-S2(y))/hbar]
      = 1-K r y^4/hbar
        + r^2[K^2 y^8/(2 hbar^2)-K y^6/(2 hbar)] + O(r^3).

Exact Gaussian moments verify invariant physical observables
F^n for n=0,2,4,6 through r^2. Omitting the Jacobian changes the
partition function at order r by -3 hbar/K; using y^2 instead
of F^2 changes its moment by -6 hbar^2/K^2. Omitting the generated
sextic changes the partition function at r^2 by
15 hbar^2/(2 K^2). These nonzero controls test the bookkeeping;
they are not numerical evidence for a continuum amplitude.

In the continuum, the physical-source identity gives the same
stable-light scattering amplitude in the full transformed theory.
The ordinary renormalized Psi has one-particle overlap
o=1-hbar w and residue Z=o^2=1-2 hbar w, where w is computed
from the actual derivative map, not set to zero. Near four
external light poles, spectral pole factorization supplies
o^4(A0+hbar A1) in the connected numerator. Amputation by the
four full two-point functions divides this by Z^4, giving

    Gamma_Psi,4 = A0 + hbar(A1+4w A0).

Proper LSZ multiplies by Z^2 and returns A0+hbar A1. The native
algebra checks all intermediate factors, not just the last
cancellation. No unstable heavy particle is used as an LSZ state.

Off-shell ordinary effective actions need not transform as scalars.
The Hessian chain rule contains S'_a F''_a as well as F'^T S'' F'.
An explicit Gaussian-plus-quartic one-variable example gives a
nonzero quartic difference -24 r^2-lambda4 r/K between the
one-loop transformed effective action (including its Jacobian)
and Gamma1(F). This is another negative control against extending
the on-shell conclusion to arbitrary off-shell 1PI functions.
