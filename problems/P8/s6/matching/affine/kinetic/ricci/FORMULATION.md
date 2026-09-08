# S6.40: antisymmetric Ricci-difference deformation

Status: exact quadratic reduction and principal-sign cases certified.
Original P8 remains OPEN.
This is a new quadratic rolling-action calculation; ancestors are frozen.

## Literal action and conventions

Keep the unrestricted S6.37 independent-connection CD/M1 action and its
original physical metric, clock and free chi. Use source signature
(-+++), derivative index last, and

    R^a_bmu_nu = partial_mu Gamma^a_bnu - partial_nu Gamma^a_bmu
                  + Gamma^a_cmu Gamma^c_bnu - Gamma^a_cnu Gamma^c_bmu,
    F13_mu_nu = R^a_nu_a_mu,
    F14_mu_nu = g_nu_r g^as R^r_s_a_mu,
    C_mu_nu = (F13_mu_nu-F14_mu_nu)_[mu,nu].

Antisymmetrization has weight 1/2. Add only

    Delta S = -(lambda_physical/4) integral sqrt(-g) C_mu_nu C^mu_nu,

with constant real coupling. This C is projectively invariant and
vanishes on every Levi-Civita connection. It is not a mixture of the
two distortion-trace curls: an explicit distortion with both traces
zero still has nonzero C. No torsion restriction is imposed.

After extracting M² tau², lambda=lambda_physical/(M² tau²),
C_normalized=tau² C_physical and q=(tau*k_physical)². The old rolling
solution has zero distortion and is preserved exactly.

## Required distinction between flat and rolling elimination

At quadratic order on that trajectory the connection Hessian is the
pure p=1/2 Palatini Hessian, not a kinematically restricted connection.
The full 64-equation contraction gives a null flat differential Schur
symbol. The corresponding curved covariant commutator is not zero:

    C[nabla Z] = (R_background/3) Z
              = 8(1+7u²)/(1+u²)² * Z

for an antisymmetric two-form Z on this conformally flat background.
The complete rolling stationary source, including coefficient jets, is

    Cstar_0i = -10u/(1+u²)^4 * partial_i n, Cstar_ij=0.

Every connection equation must be retained before eliminating Z on
1+lambda*R_background/3 != 0. No flat-space pole or raw lapse derivative
is substituted for the curved constrained action.

## Scientific scope

Derive the exact reduced quadratic addition a³*q*A(u,lambda)*n²,
including its sign, full scalar constraints and time-dependent spatial
boundary. Test positive principal kinetic and gradient matrices for
lambda>=0; a gyroscopic term is retained. Treat negative couplings,
singular algebraic charts, the center Theta=0 chart and lambda=0
separately. The regular first-order center action, not a division by
Theta, supplies the crossing interface.

A positive-principal quadratic deformation is not an all-frequency or
nonlinear stability theorem. A nonzero change to the physical CD/M1
quadratic action is not exact matching. A momentum-limited relative
kinetic bound is not a justified EFT frequency cutoff or a full V/G/B
remainder. The calculation supplies no additional heavy spectrum,
vacuum or finite-gravity estimate, UV completion, original P8 closure,
or new requirement for the completed scoped photon objective.
