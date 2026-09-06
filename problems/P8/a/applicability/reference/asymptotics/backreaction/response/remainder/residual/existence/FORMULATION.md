# A.10: actual nonlinear response continuity and a causal inverse

This gate supplies two quantitative analytic ingredients for a possible
actual semiclassical solution theorem. It does **not** prove that theorem.
The pinned A.9/A.8 metric, state and finite prescription remain unchanged.

## Actual mode functional

For real `C1` potentials `U,V` on `[0,T]`, with `U(0)=V(0)=0`, use the
same original plane-wave data for the exact modes
`w_k''+(k^2+U)w_k=0`, `w_k(0)=1`, `w_k'(0)=-ik`. Let
`F_k=|w_k|^2=2k|v_k|^2`. Subtract only the exact degree-zero and
degree-one terms of its potential expansion and define

    R[U]=integral_0^infinity k*(F_k[U]-1-F1_k[U]) dk.

For smooth potentials flat in the preparation past this is the actual
nonlinear mode part of A.7's renormalized Wick square, with coefficient
`hbar/(4*pi^2*a^2)`. The quadratic integral is defined by its explicit
Abel limit; the higher-order integrals converge absolutely. The same
formulas define a `C1`-valued continuous nonlinear functional on the stated
`C1` domain. They do not assign a smooth Hadamard metric or quantum state
to every element of that domain, or infer uniqueness of extension merely
from the smaller smooth past-flat subset.

Put `D=||U'-V'||`, assume `||U'||,||V'||<=M`, and set `z=M*T^3`.
Then

\[
 \|R'[U]-R'[V]\|\le(2z+18z^2e^{2z})D,
 \qquad \|R[U]-R[V]\|\le T(2z+18z^2e^{2z})D.
\]

If the potentials agree through `T-L`, `0<L<=T`, the derivative
coefficient improves to

\[
 C_{\rm shared}=MTL^2\left[\frac54+\frac12\log(T/L)\right]
                   +18M^2T^5L e^{2MT^3},
\]

and the value coefficient is `L*C_shared`. This retains the full common
nonzero prehistory. No reset to a different state at the future starting
slice, positive frequency gap, or higher potential derivative is used.
The [proof](notes/mode-response.md) includes all Dyson insertion,
free-phase cancellation, initial endpoint and quadratic Abel arguments.

The analytic exponential estimates hold for finite `z`. The numeric API
uses the explicit subdomain `0<=z<=1/4` and rational exponential/logarithm
upper bounds. At dimensionless `delta=10^-14`, history bound `T=3`,
`M=2*delta*(61013499/8192)` and future length `L=1/4`, the derivative
coefficients are respectively below `10^-8` and `10^-10`. A.7 supplies
the actual derivative bound used to place its prepared potential in this
ball. These are coefficients for `R` alone, not complete SEE map norms.

## Isolated retarded logarithmic inverse

For real `beta`, the dimensionless causal operator is

\[
 D_\beta f=\frac1{8\pi^2}
  \left[(\beta-\gamma_E)(f-f(0))
       -\int_0^t\log(t-s)f'(s)\,ds\right].
\]

Its inverse kernel is

\[
 K_\beta(t)=8\pi^2\left[e^{-\beta}e^{e^{-\beta}t}
  +\int_0^\infty\frac{e^{-rt}}{(\log r+\beta)^2+\pi^2}\,dr\right]>0.
\]

The [inverse proof](notes/inverse.md) keeps both the positive pole and
the branch cut, establishes the exact `C0` convolution norm and a finite
bound tending to zero with interval length. A rational dyadic example
has norm `<1/3` for `beta>=-1` and `L=2^-1024`. It is only an operator
calibration: it does not designate a physical proper-time interval.
`beta` is not an independently selected finite-curvature counterterm.

An explicitly conditional Banach check requires a separately established
full nonlinear Lipschitz constant, fixed-point residual, positive ball
radius, `q<1`, and the self-map inequality. The illustrative rational
inputs are not computed SEE data. Neither this inverse nor the mode
lemma proves a `C1`-to-`C1` inverse or automatic smoothness at the initial
surface.

## Evidence boundary

The [certificate](certificates/response-inverse.json) pins A.9 and replays
its lineage; independent Fraction arithmetic reconstructs the constants.
Independent phase-product, radial, endpoint and source-normalization
audits are included. The analytic kernel, convergence and functional
arguments are written proofs, not Lean formalization.

The [remaining SEE obligations](notes/see-gap.md) are explicit: the full
positive-metric map and its local geometric terms, exact compatible
initial constraint and joint state preparation, a complete interval
enclosure, and a smoothness/Hadamard argument are not supplied. No nearby
exact solution, perturbation of A.9 to an unknown solution, massive or
interacting field extension, cosmological incompleteness, or full P8
closure follows.
