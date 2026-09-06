# P8(a) A.9 — Actual perturbed all-sampler QSEI

Adopted 2026-09-06. The immediate immutable dependency is A.8, certificate
SHA256 `dbd88193542224b3eeb76e18923744f3317791ad46961caca3562a1525ba56d1`.
No earlier metric, state, prescription, theorem or certificate is changed.

## The theorem

Keep precisely A.8's smoothly prepared metric, one free real massless
minimally coupled scalar, and the reference Hadamard state transported
from its unchanged radiation past. Keep its numerical raw-H_lambda
prescription `lambda=2*sqrt(2)*A*eta_star²`, physical gamma=0, no additional
finite Wick-square shift, and

    0<=delta=16*epsilon*d/(A²*eta_star⁴)<=10^-14,
    d=kappa*hbar/(46080*pi²), y=sqrt(2*t/A)/eta_star.

For EVERY Hadamard target state omega of this same field algebra and EVERY real smooth proper-time
sampler h compactly supported in the target 2<y<3,

    integral h² * [<T_uu>_omega-<T>_omega/2] dt
        >= -5*hbar/(16*pi²)*integral |h_ddot|² dt.

Here u is the unit comoving timelike velocity in the FK +--- convention.
The target state need not be homogeneous, isotropic, quasifree, zero-mean
or have bounded pointwise energy. The statement extends by density to
real H2_0 samplers on each compact interval inside the target. Spatial
homogeneity of the reference geometry and state gives the same bound
along every comoving line. It does not assert coverage of noncomoving
geodesics or of times outside the target envelope.

The coefficient 5 is a convenient sufficient constant, not an optimum.
The sharper derived bound is `(2026/961+delta*Rroot)²`, where Rroot is
the exact rational reported by the certificate and Rroot<64466343.
It is below 5 throughout the closed amplitude interval. The proper-time
coefficient contains no arbitrary normalization A or conformal length;
physical epsilon=1 is included in A.8's explicit preparation-time regime.

## Proof obligations discharged

The bound is derived from the actual pulled-back two-point function, not
inferred from the coincident stress remainder. Phase-resolved exact modes
have uniform derivative bounds and inverse-frequency bounds through
three derivatives. Flat past jets are essential to the retarded derivative
transfer. The error in the differentiated modes has two uniform derivative
bounds, sufficient for two sampler integrations by parts and an integrable
joint ultraviolet k/alpha moment. The infrared integral uses full complex
Parseval and a nonsingular Volterra estimate.

A separately computed auxiliary flat kernel retains the actual conformal
Hubble coefficient; it is only a spectral norm-comparison device, not an
instantaneous vacuum or positive perturbed-metric state. Exact proper-time
conversion and Poincare inequalities give a bound uniform over all samplers
in the stated domain. Finally, A.8's proved positive actual reference EED
turns the renormalization-independent difference bound into the displayed
absolute bound in the specified scheme.

The certificate replays A.8 and its entire pinned ancestry, exact phase,
spectral and clock identities, positive rational bounds, a separately
assembled Fraction-only jet-polynomial recurrence, and adversarial tests.
The positive-type/Volterra/Sobolev arguments remain written proofs, not
Lean FORMALIZED theorems. Finite regressions alone do not prove the
all-mode or all-sampler quantifiers.

## What remains open

The prepared metric is off shell. Neither A.8's small SEE residual nor
this actual-metric field inequality supplies an exact or nearby
self-consistent SEE solution. The whole geodesic domain and initial
curvature conditions of A.1 have not been discharged. Massive,
nonminimal and interacting realistic fields, arbitrary finite schemes,
other metric families, and cosmological-scale incompleteness remain
outside the theorem. This closes the perturbed all-sampler input for
this fixed testbed, not P8(a) or P8 overall.
