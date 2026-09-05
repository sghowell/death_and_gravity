# P8-S5.6.CD: fully spatially reduced M1 interactions

This checkpoint advances the fixed CD witness with rolling free canonical matter
from S5.5 invariant lapse reduction to physical Fourier interactions. It leaves
all published input files and certificates unchanged. It does **not** complete
P8(b), normalize the coupled propagating modes, or establish interaction control.

## Frozen scope

Use exactly P8-2.CD and S5.5: the same covariant tube, signature, physical matter
metric, background, nonlinear auxiliary spatial metric and local canonical
primitive. In dimensionless cosmic time `u`, `d=1+u^2`, the background has
`a=d^2` and `chi_dot=1/(10*d^6)`. The extra scalar is free and minimally coupled
to the original physical metric, not redefined to another causal frame.

The spatial gauge is the small-field formal chart
`hat_h_ij=a^2*g_ij`, `g_ij=(1+2*zeta)*delta_ij+gammaTT_ij`.
Work through perturbative degree four. Retain both transverse-traceless tensor
polarizations, the gravitational scalar and the matter scalar, and all three
matter-sourced spatial momentum constraints. This is a formal local reduction;
no global nonlinear gauge, zero-mode backreaction, Dirac or PDE theorem is added.

External momenta are arbitrary rational vectors in three dimensions. Their total
is zero, and every nonempty proper subset must have nonzero total momentum.
TT polarizations are arbitrary nonzero real rational symmetric matrices obeying
`tr E=0` and `E*k=0`. Distinct labelled nilpotent amplitudes extract the coefficient
of their product, with no inserted `1/n!`; this includes repeated field species.

At evaluation use fixed local units `ell=tau*sqrt(d)` and spatial `a=1`.
Let `q=ell^2*k_physical^2` and compactify time by

`r=u/(1+sqrt(1+u^2)), x=2*r/(1+r^2), y=(1-r^2)/(1+r^2)`.

Finite physical time has `-1<r<1`; `r=+-1` denotes the two limiting tails only.
The background functions in this note are the scaled functions

`H=4*x`, `l=y^11/10`, `Theta=x*(4-y^6)`, `Lambda=1-3*y^6/2`,

`w=l*(3*y^6/2-1)`, `J=ell^2*J_physical`, `J0=J+w^2/2`.

The pinned S5.5 proof gives `1/10<J<8` on finite time and the compact extension.

## Claim

In this scope, the finite York recursion solves the three spatial constraints
through degree three. Substituting its result and full metric/matter geometry
into the pinned stationary lapse formula gives the cubic and quartic canonical
phase kernels. The required mixed scalar-matter boundary and all background
time-generator terms are retained. The same construction has a regular phase
gamma chart, without inverse `Theta` or `Lambda`.

For fixed nonexceptional external momenta in these local units, these phase coefficients extend
continuously to compact time and are bounded there by the algebraic closure
argument in [notes/reduction.md](notes/reduction.md). This is not an explicit
uniform momentum bound, a velocity bound, or an interacting cutoff estimate.
Higher-order report fixtures use exact rational times; no expanded symbolic
quartic majorant is claimed. The full quadratic bridge is checked at symbolic
time, in both charts, against an independently derived action.
This is not a uniform statement along a fixed-comoving-mode history, whose
local q tends to zero in the two tails.

In admissible velocity charts, the full two-scalar quadratic inverse Hessian
gives the **unnormalized** velocity substitution and the quartic Legendre contact
`H3_P^T*A^-1*H3_P/2`, including scalar off-diagonal terms and both tensor channels.
This is not the propagating cubic-exchange contribution to a scattering amplitude.

## Velocity domains and the finite-momentum distinction

The ordinary scalar chart requires `Theta!=0`. The positive gamma velocity chart
requires `q*Lambda^2>J0`, at every external and internal momentum used in the
quartic contact. A symbolic-time API call only returns a rational expression;
it does not prove that its velocity chart covers the whole compact interval.
Numeric calls reject singular or non-positive kinetic matrices.

At the bounce `J=1199/800`, `J0=3/2`, `Lambda=-1/2`, `w=1/20`.
The finite-momentum gamma kinetic matrix at `q=8` is

`K=[[24,1/5],[1/5,401/800]]`.

Its high-frequency limit is the earlier principal matrix
`[[6,1/20],[1/20,1/2]]`. In particular, the finite-q matter diagonal is generally
not `1/2`. The S5.5 principal-cone theorem must not be applied to that finite-q
matrix. The gamma velocity chart is singular at bounce `q=6`, while the phase
Hamiltonian is finite; this alone is not a physical instability claim.

## Remaining work

The next CD/M1 step is a coupled, time-dependent quadratic mode normalization
and free-evolution estimate with all induced derivative/connection terms. In
moving local units the fixed-comoving mode obeys `Dq=-6*x*q`, unlike D-only M0.
These terms must then enter physical interaction majorants and a specified
finite-time hard-transfer tree criterion. None is inherited merely from S5.3.D.

No canonical configuration/mode normalization, finite-k oscillator positivity,
interaction window, inclusive amplitude, all-orders/loop estimate, radiative
protection of the exceptional relation, nonlinear/BKL stability, or UV
admissibility/completion is certified here.
