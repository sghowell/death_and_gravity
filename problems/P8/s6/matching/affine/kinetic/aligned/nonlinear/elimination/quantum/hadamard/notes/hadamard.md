# Proca Hadamard comparison for the new all-order Cauchy state

The state here is the separate cutoff-summed preparation of
[construction.md](construction.md). The old finite-order state is
not relabelled. The statement concerns the retained free Gaussian
vector on the actual clock, where its mass is constant and its
field equation is ordinary massive Proca. It does not assert an
interacting UV completion or a state for every field in that model.

## A constrained-vector Hadamard reference

Let chi be the smooth turn-on used in the construction and set

    a_aux(u)=1+chi(u+3)[a_actual(u)-1].

This auxiliary scale factor is one for u<=-2 and equals the actual
(1+u^2)^2 for u>=-1. It is smooth and at least one. Both its
FLRW spacetime and the actual one are globally hyperbolic with
Cauchy surfaces R^3: causal coordinate speeds on a finite time
interval are bounded, and an inextendible causal curve cannot
escape spatial infinity in finite u. This auxiliary geometry is
used only to construct a reference state. It is not a replacement
for the physical bounce or a solution of its light equations.

In the flat region use the positive-frequency Minkowski Proca
vacuum. The code independently checks its three positive spacelike
polarizations, mass shell, Lorenz constraints and polarization sum.
For momentum along the third axis the two transverse polarizations
are the first two unit spatial vectors and the longitudinal
covector is (-k/m,0,0,omega/m). Their Gram matrix is the identity
and their sum is g_munu+p_mu*p_nu/m^2. This is a positive
three-mode Gaussian state on the constrained Proca algebra, not
four unconstrained scalar states with a positive fiber metric.
Its two-point distribution has the usual positive-frequency
massive scalar singularity with the Proca differential projector,
hence the Hadamard wavefront set.

The Proca time-slice and Hadamard propagation property extends
this reference from a flat Cauchy neighborhood near u=-3 to
the auxiliary spacetime. Restrict near u0=-1/2, where the metrics
coincide, then extend to the actual spacetime. This uses
[Moretti, Murro and Volpe, Proposition 4.7](https://arxiv.org/html/2210.09278v3).
Its hypotheses are the globally hyperbolic Proca spacetimes and
Cauchy neighborhoods just specified, with fixed m^2>0.
Spatial translations and rotations are preserved by both
geometries and the initial vacuum, so the reference has the
same transverse/longitudinal Fourier decomposition.

## All-order comparison with that reference

On the compact auxiliary interval [-2,u0], repeat the oscillator
recurrence with its smooth Hubble function. For each fixed N all
coefficients and their finitely many derivatives are bounded.
At sufficiently large k, omega is uniformly comparable to k,
the finite W_N is positive, and its reference residual is
O(k^(-2N)). These are qualitative estimates for this comparison,
not constants used in the physical error budget.

All Hubble jets vanish at the flat endpoint u=-2. Thus the W_N
mode there has exactly the Minkowski initial data, for every N.
The normalized reference-basis variation equation on the compact
interval then gives a mode-mixing error O(k^(-2N-1)) relative to
the auxiliary exact reference. The coefficient norm estimate
uses residual divided by frequency; it does not bound the
oscillator by an exponentially growing exp(constant*k) norm.
This is the same exact variation-of-constants identity as in
S6.50, applied at arbitrary fixed N.

At u0 the new cutoff-summed frequency and its slope agree with
W_N,W_N' up to O(k^(-2N-1)), by the all-order cutoff estimate.
The exact Cauchy Bogoliubov formula therefore gives error
O(k^(-2N-2)) relative to that W_N preparation. Composing the
two comparisons shows that mixing between the new state and
the propagated Proca Hadamard reference decreases faster than
every inverse power of k. Arbitrary overall phases do not
affect the Gaussian two-point function.

Both states evolve by the same exact physical Proca equation
after their common Cauchy slice. On any compact actual time
interval, canonical mode values and any fixed finite number of
time derivatives have at most polynomial growth in k. For
example, a positive-frequency oscillator energy estimate at
large k has a k-independent finite-interval exponent because
the coefficients and their normalized logarithmic derivatives
are bounded there. Higher time derivatives then follow directly
from the smooth second-order equation. The rapid Cauchy mixing
therefore remains rapid enough after every fixed number of
time derivatives.

Reconstruct the physical vector, including its temporal
constraint. The transverse multiplier is 1/sqrt(a). For a
longitudinal canonical mode, the spatial multiplier has magnitude
omega/(sqrt(a)*m), and the temporal component is

    k*(v'-d_L*v)/(a^(5/2)*m*omega)

up to the Fourier phase convention. These identities are checked
against the original constrained canonical normalization. The
multipliers and all fixed time derivatives grow at most
polynomially in k on compact time intervals. Their derivatives
cannot destroy a faster-than-every-power comparison.

Thus the Fourier integral for the difference of the two physical
vector two-point functions converges absolutely after any fixed
number of spacetime derivatives, uniformly on compact sets.
The bounded low-momentum part causes no difficulty: m>0, and
the physical constrained readouts are regular; near k=0 the new
state is exactly the frozen state because the cutoffs vanish.
The difference is a smooth bidistribution.

The new positive normalized three-polarization Gaussian therefore
has the same Hadamard wavefront set as the reference. This proves
the all-order state property for the actual retained Proca field.
It does not rely on importing a scalar low-energy-state theorem,
or on assuming that fourth-order Cauchy data already satisfy
the all-order condition.

## Quantitative versus qualitative conclusions

The auxiliary reference is used only to establish the all-order
short-distance property. No bound on its unknown energy or on
the constants in this qualitative asymptotic comparison is
inserted into the physical estimates. Those estimates instead
use the explicit C_H/nu^6 initial comparison with S6.50 and the
unchanged local prescription, as derived in the transfer proof.

Hadamard regularity here is a property of the selected free
Gaussian component. Quantitative higher-time-derivative bounds,
quantum functional variations, other fields and interacting
loops, corrected light equations and cones, cutoff, vacuum/Regge
matching and the adopted V/G/B conditions remain separate work.
