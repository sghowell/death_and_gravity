# Literal Ricci difference and the full curved connection solve

## 1. Source and unrestricted convention

The Ricci contractions are the ones in equations (2), (3) and (12) of
[Barker and Marzo, arXiv:2402.07641v2](https://arxiv.org/html/2402.07641v2),
translated explicitly to the derivative-last source convention in the
formulation. Only that definition motivates this new action. No
Minkowski spectral-health result, torsion-free assumption or proposed
complete classification is imported from the paper.

With kappa=Gamma-LC(g), put

    N[kappa]_mu_nu = [nabla_a(kappa^a_nu_mu
                         -g_nu_r g^ab kappa^r_b_mu)
                       -nabla_mu(V_nu-U_nu)]_[mu,nu].

This is the linearized C tensor at zero background distortion.
The Levi-Civita part of C vanishes identically for arbitrary g, so its
metric variation does not supply an omitted independent term here.
At quadratic order the curvature-product terms quadratic in kappa are
irrelevant to C², since the background C and kappa both vanish.

The full six-by-64 flat symbol N(k) is constructed from these literal
curvature contractions with all four covector components independent.
It annihilates all four projective columns. A dense unrestricted
curvature calculation and an exact Lorentz transformation check the
contraction independently. For k=(1,0,0,0), kappa^0_12=1 and
kappa^0_21=-1 give C12=-1 while V=U=0. Thus this operator is outside
the already closed constant two-trace family.

## 2. All sixty quotient equations and the null flat symbol

Let M be the frozen 60-dimensional projective-quotient Hessian, E its
embedding in the unrestricted 64 components, and Nq=N E. On the full
actual rolling trajectory p=1/2, the coefficients c and J4 in the
connection-quadratic action vanish. M is therefore the pure Palatini
quadratic form. Coefficient variations still contribute to the source
and are not set to zero.

All original quotient blocks are inverted directly. With W2 the
six-component two-form pairing diag(-1,-1,-1,1,1,1), the exact result is

    Nq M^-1 Nq^T = 0,
    [E M^-1 Nq^T W2 Z]^a_bc = k_c Z^a_b.

The first equality holds coefficientwise in every component of k; it
is not just a homogeneous limit. The second is checked against all
64 Euler components as well as the quotient projection. Nq has rank
six on timelike, spacelike and null nonzero covectors. Lorentz
covariance and scaling cover all three orbits. The flat matrix update
is nilpotent after multiplication by M^-1, yielding an exact polynomial
inverse and unchanged determinant. It creates no new connection pole
in this constant-coefficient background. This fact alone does not
settle the rolling metric/clock dynamics.

## 3. The curved commutator must not be discarded

The tensor identity for the lift extends covariantly because the
on-clock Hessian depends only on the background metric, whose
Levi-Civita derivative is zero. In covariant variational notation,

    M nabla Z = -N^dagger W2 Z.

The minus sign comes from the adjoint of a derivative, with the full
background volume element included. For kappa^a_bc=nabla_c Z^a_b,
where Z_ab is antisymmetric, V=-U=nabla_a Z^a_mu. Hence

    C[nabla Z]_mu_nu
      = [nabla_a,nabla_mu] Z^a_nu
          -[nabla_a,nabla_nu] Z^a_mu.

Evaluating the commutator gives

    Ric_r_mu Z^r_nu - R^r_nu_a_mu Z^a_r - (mu <-> nu).

For flat FLRW, R^0_i0j=(H_dot+H²)delta_ij in an orthonormal frame
and R^i_jkl=H²(delta_ik delta_jl-delta_il delta_jk). Direct contraction
on all six two-form components gives the scalar endomorphism

    K_curv = 2(H_dot+2H²) = R_background/3.

An independent coordinate calculation starts with six arbitrary
functions Z_mu_nu(u,x,y,z), computes every component of nabla Z with
the actual FLRW Christoffels, and then applies the literal Ricci
contractions. It reproduces K_curv Z with no surviving derivatives
of Z. The actual a=(1+u²)² gives

    K_curv = 8(1+7u²)/(1+u²)² > 0

for every finite real time. The flat cancellation was commuting
derivatives; dropping the curved remainder would give a wrong answer.

## 4. Exact quadratic multiplier reduction

Let y=kappa-kappa_star, where the full original stationary connection
is used. Introduce an antisymmetric two-form Z before eliminating
connection components. In pair notation, the new quadratic density is

    L = L_CD/M1 + y^T M y/2
          + Z^T W2 Z/(2 lambda) - Z^T W2(Cstar+N y).

This representation is used only for lambda!=0; lambda=0 is the
original auxiliary action. The y equation and the preceding full
adjoint identity give y=-nabla Z, retaining every connection equation.
Covariant integration by parts gives

    integral (nabla Z)^T M(nabla Z) = -integral Z^T W2 K_curv Z.

Consequently the complete reduced multiplier action is

    L = L_CD/M1 +(lambda^-1+K_curv) Z^T W2 Z/2-Z^T W2 Cstar.

There is no derivative of Z left. On 1+lambda K_curv!=0,

    Z=lambda Cstar/(1+lambda K_curv),
    Delta L=-lambda Cstar^T W2 Cstar/[2(1+lambda K_curv)].

The algebraic Euler equation and the action substitution are checked
independently with six unconstrained two-form components. For positive
lambda the denominator is at least one on the entire trajectory.
No homogeneous connection state or retarded inversion remains on
this quadratic branch. The Lorentzian pairing is algebraic and is
not by itself a propagating-kinetic sign test.

These are quadratic statements at the actual solution, not a nonlinear
open-tube elimination: higher-order curvature/distortion couplings and
other backgrounds require their own rank and field-equation analysis.

## 5. Explicit off-clock nonclaim control

The generic original CD quotient can also be contracted before putting
p=1/2. For k=(omega,0,0,k_z), the exact constant-coefficient result is

    Nq M(p)^-1 Nq^T = k_z² diag(-a_p,-a_p,0,0,b_p,b_p),
    a_p=(2p-1)²(3p+1)/[2p(8p²-1)], b_p=(2p-1)²/(8p²).

Both are positive on the original closed p tube away from p=1/2.
The response has rank four for k_z!=0 there, yet contains no omega.
Two independent exact p values inside that tube verify this control.
It proves that the null identity used on the actual clock is not a
nonlinear open-tube identity. Its purely spatial character also means
that a nonzero Schur response alone is not evidence of a new time mode
or a ghost. General variable-coefficient inversion and the complete
metric constraints remain necessary for any such separate conclusion.
