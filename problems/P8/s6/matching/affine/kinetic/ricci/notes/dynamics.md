# Actual source, constrained principal matrices, and matching limit

## 1. All stationary connection jets before contraction

The input is the original full S6.37 stationary connection, not an
isolated connection ansatz. On the clock, p=1/2, x=-1, h=(1+u²)^3,
and h_dot/h=3H/2. The lower-coefficient ODE supplies q=q_X=0 there.
In the timelike orthonormal scalar frame, the complete lapse jets are

    delta p=-n/(2h), delta s_clock=-n,
    delta p_phi=delta J3=h_dot*n/(2h²),
    delta p_X=-n/(4h²), delta c_phi=h_dot*n/h²,
    delta c_X=2n*c_XX, c_XX=-(4h-1)/(4h²).

They are applied to every component of the generic frozen solution.
The orthonormal Hessian jets are H00=-n_dot, H0i=-D_i n, while all
six spatial Hessian perturbations are initially independent. The
spatial shift contribution to the mixed Hessian cancels only after
using the normal frame, as in the full frozen ADM calculation. The
six spatial perturbations then drop out of the complete connection
variation, leaving the following nonzero components:

    kappa^0_00=-15Hn/(8h),
    kappa^0_i0=kappa^i_00=-3D_i n/(2h),
    kappa^0_ij=kappa^i_0j=delta_ij*(Hn-n_dot)/(2h),
    kappa^i_j0=delta_ij*5Hn/(8h),
    kappa^i_jk=(delta_jk D_i n-delta_ik D_j n)/(2h).

All 64 equations are checked against this expression, including its
zero entries. The exact background stationary connection is zero.
Therefore perturbing the coframe multiplies zero, and conversion to
coordinate components uses only the background coframe. The literal
curved Ricci contraction, with every FLRW Christoffel retained, gives

    Cstar_0i=r(u)*partial_i n, Cstar_ij=0,
    r(u)=-5H/(2h)=-10u/(1+u²)^4.

The lapse-time derivatives cancel. Keeping only a flat principal
symbol would miss this nonzero rolling source; keeping only selected
second lapse derivatives would falsely introduce an acceleration.

## 2. Actual reduced quadratic correction

The complete curved multiplier solve in the geometry note yields

    Delta L/a³ = q*A*n²,
    A=lambda*r²/[2(1+lambda*K_curv)]
      =50lambda*u²/[(1+u²)^6*((1+u²)^2+8lambda*(1+7u²))].

The sign follows from Cstar having only electric components in source
signature (-+++), with Cstar_mu_nu Cstar^mu_nu=-2q*r²*n².
For lambda>=0 the denominator is strictly positive for every time,
and A>=0. At u=0, A=A_dot=0. The new terms preserve the background
and supply no independent tensor, metric-vector or connection modes
on this regular quadratic branch: the full source is scalar and the
six two-form variables have been eliminated algebraically.

## 3. Physical configuration and crossing constraints

Add q*A*n² to the complete old metric/free-chi scalar action, retaining
lapse n and shift b before any elimination. The shift equation remains

    v_dot=Theta*n-ell*s/2, s=delta chi.

For finite u!=0, Theta!=0, the lapse is fixed by this equation and its
companion equation reconstructs b. The physical velocity form is

    L_kin/a³ = (s_dot+w*v_dot/Theta)²/2
                  +(J+q*A)*v_dot²/Theta².

The determinant is (J+q*A)/(2Theta²). The original J has an even
nonnegative numerator with positive constant and positive denominator.
Thus both scalar kinetic directions are positive for lambda>=0 and
q>0 throughout this punctured chart. Independent joint solves of the
literal n,b equations reproduce this matrix at both time signs and
nonunit coupling; no metric or matter fluctuation is frozen.

The center is handled instead with the regular first-order density

    H=C0+R0²/(J+q*A),
    R0=q*(Theta*b+Lambda*v)+w*pm/2-3ell*s*Theta/2,
    C0=pm²/2+q*ell*b*s+q*s²/2-q*v²-3ell²*s²/4.

It is derived by the full Legendre/shift reduction before eliminating
the lapse, whose solution is -R0/(J+q*A). This expression has no
inverse Theta. The canonical exchange P_b=2a³qv uses the same moving
boundary -Hubble*b*P_b as the original theory. For lambda>=0 all its
coefficients are smooth at the crossing. Since A=A_dot=0 there, its
Hamiltonian and first coefficient jet agree with the old first-order
system. The old regular center kinetic chart remains a valid control
at q>6. No u=0 substitution in the divergent v-velocity chart is used.
The direct canonical momentum Hessian has determinant
100(q-6)/(1199q), positive matter pivot 1200/1199, and inverse kinetic
matrix

    [[6q/(q-6), q/(20(q-6))],
     [q/(20(q-6)), (200q-1199)/(400(q-6))]].

These identities are replayed from the new first-order action with its
actual coefficients and the moving canonical boundary.

## 4. Time-dependent spatial boundary and gyroscopic mixing

Kinetic positivity alone is insufficient. Let

    B=A*ell/Theta², a=(1+u²)²,
    Gvv=(1/a)*d_u(a*Lambda/Theta)-1,
    Gss=1/2-A*ell²/(4Theta²),
    Gvs=-Lambda*ell/(2Theta)+(1/(4a))*d_u(a*B).

The exact q-dependent scalar action, after its time boundary, contains

    q*A*v_dot²/Theta² + q*B*(s*v_dot-v*s_dot)/2
        -q*(Gvv*v²+2Gvs*v*s+Gss*s²).

The middle term is retained, not removed by pretending B is constant.
The boundary primitive divided by the constant comoving momentum
squared is a*[(Lambda/Theta)*v²+B*v*s/2]. The code differentiates it
with the actual time dependence and checks the entire action identity.
In particular a³q is proportional to a, not to a³, since q=k_com²/a².

The symmetric matrix G is positive definite on every finite u!=0,
lambda>=0. This is a continuous exact polynomial proof, not a scan:

- Gss has a numerator with 37 monomials, degrees (36,1) in (u,lambda),
  positive constant 36, only even u powers and nonnegative coefficients.
- det(G) has 180 numerator monomials, degrees (74,4), positive constant
  97119, only even u powers and nonnegative coefficients.
- Their denominators are positive on the stated domain, with the
  punctured factor u² separated explicitly in the determinant proof.

Gvv itself has negative polynomial coefficients; no incorrect
all-positive-coefficients argument is made for that entry. The strict
Gss and determinant tests establish the two-dimensional inertia.

For the principal action L=velocity^T K velocity+velocity^T Q*x
-q*x^T G*x, Q is antisymmetric and includes the displayed gyroscopic
term. Its Hamiltonian is

    H_principal=(p-Q*x)^T K^-1(p-Q*x)/4+q*x^T G*x.

An independent generic Legendre check verifies this square. Positive
K and G therefore survive the gyroscopic mixing. This is a principal
quadratic kinetic/gradient statement, not conservation of an energy
with time-dependent coefficients, an all-frequency mode bound or
nonlinear stability. The original tensor sector is unchanged.

## 5. Negative coupling and singular-chart controls

Every lambda<0 has a regular tail interval with 1+lambda*K_curv>1/2:
use u²>112*abs(lambda), since K_curv<56/u². There H is nonzero and
A<0. Choosing q>J/(-A) gives one negative physical scalar kinetic
direction after both constraints. This rejects global literal
kinetic health for every negative coupling without treating an
algebraic pole as an automatic physical instability.

Other times can have 1+lambda*K_curv=0. Such charts are not inverted.
Exact tests include lambda=-3/49 at u²=5/7 and lambda=-1/8 at u=0
or u²=5. No continuation through those rank changes is asserted.
Lambda=0 is the unchanged auxiliary action, handled separately from
the multiplier representation that divides by lambda.

## 6. Nonzero matching defect and its limited bound

For lambda nonzero on a regular time interval with H!=0, the physical
quadratic action changes by q*A*n². Its lapse Hessian is 2q*A, not a
boundary term. After the full physical constraints its velocity
matrix changes by q*A/Theta² in the v direction. This is not exact
CD/M1 matching in the unchanged metric/clock/matter frame.

One uniform comparison is available for lambda>=0. The polynomial
numerator of 10-r²/(2J) is strictly positive on all real time, so

    0 <= q*A/J < 10lambda*q    (lambda>0,q>0).

In the original completed kinetic squares this is the only nonzero
relative eigenvalue of the added kinetic form. On a declared momentum
range q<=Q, 10lambda*Q<=epsilon gives
K_original <= K_new <= (1+epsilon)*K_original on the punctured chart.
At lambda=0 the two matrices are identical. In physical units the
bound is 10lambda_physical*k_physical²/M², with tau cancelling.

This is only a quadratic kinetic-coefficient bound. It is neither a
full action/gradient/loop remainder nor a heavy-state or frequency
cutoff estimate. M²=3,tau=2,lambda_physical=5 checks normalized
lambda=5/12, C_normalized=4C_physical, and at u=1/2 the curved
denominator is 103/15. Dropping the curvature commutator fails that
nonunit control. No extra mass or heavy gap is assigned.

No nonlinear open-tube inverse, vacuum extension, finite-gravity
dispersion estimate, V/G/B closure, UV completion or original P8
closure follows. The old photon and linear CD/M1 certificate scopes
are unchanged.
