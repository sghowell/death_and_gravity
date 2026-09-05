# S6.2 — Metric-map matching and the global tail obstruction

This extends the [first parent tests](parent-tests.md) without changing the
old witness. It separates two statements: an Einstein/NEC parent cannot
match the required geometry and tensor tails by a regular disformal map;
and even the broader route of mapping the exact D tube to quartic Horndeski
loses invertibility. Neither is a no-go for all parent theories.

## 1. Derive the tensor matching dictionary first

For homogeneous clock fields, write the auxiliary parent metric as
`g_E=C*g+B*dphi*dphi`, C>0, W=C+BX>0. Physical time is still t and

```
a_E=sqrt(C)*a,  dT_E=sqrt(W)*dt.
```

For an Einstein seed with positive constant Planck mass M, pull back its
tensor action `M^2/8 integral dT_E a_E^3 [(dh/dT_E)^2-(grad h)^2/a_E^2]`.
Keeping all measure and derivative factors gives

```
g_T=G_T/M^2=C^(3/2)/sqrt(W),  f_T=F_T/M^2=sqrt(C*W),
C=sqrt(g_T*f_T),  W=f_T^(3/2)/sqrt(g_T).
```

The two tensors and their inverse are replayed independently. Minimal scalar
kinetic terms contain no tensor derivatives on a homogeneous background.
For a general seed, the same kinematic factors multiply its own G_T,E and
F_T,E; setting them both to M^2 is an additional Einstein-seed assumption.
Canonical tensor normalization alone says nothing about the scalar NEC.

The conformal case W=C=q recovers G_T=F_T=M^2*q. Exact D tensor matching
G_T=F_T=M^2 forces C=W=1 on the trajectory; an Einstein/NEC seed then has
the same scale factor and cannot bounce. Equality on the entire open tube
is stronger: C=W=1 as functions there, so no X-dependent metric derivatives
can generate D's nonzero lapse-velocity mixing. This exact EH-plus-P seed
subcase fails even without assuming its P-sector satisfies NEC. These
statements concern this operator/field map, not arbitrary heavy elimination.

## 2. A global argument that does not assume small frame derivatives

Let y(T_E)=log(a_E). A classical flat Einstein solution with positive
scalar field metric satisfies

```
y''(T_E)=dH_E/dT_E=-(rho_E+p_E)/(2M^2)<=0.
```

Thus y is concave on every regular proper-time interval. Pick any finite
physical time t0. If a_E grows without bound in both physical-time tails,
there are finite t_minus<t0<t_plus with a_E(t_minus),a_E(t_plus)>a_E(t0).
Since sqrt(W)>0, proper-time order is preserved and both intervals have
positive finite length. The mean-value theorem gives a negative H_E at
some earlier time and a positive H_E at a later one, contradicting its
nonincrease. Equivalently, concavity places the middle y above its weighted
endpoint chord, which is impossible when both endpoints exceed it. The
replay checks the equivalent secant algebra; this all-time calculus argument
is a written lemma, not a conclusion from a sample grid.

For P8, a grows in both tails while positive constant tensor limits imply
`C=sqrt(G_T*F_T)/M^2` tends to a positive constant on each end. Consequently
a_E grows in both tails, and this whole Einstein/NEC parent route is
excluded. Arbitrarily large but smooth changes of C or W near the bounce
do not evade it. No assumption that C_dot or C_ddot tend to zero, and no
separate assumption about infinite Einstein proper time, is needed.
Nonminimal Q(Phi)R parents with healthy Einstein fields are included when
their actual tensor coefficients have the required tails. A lower bound
on Q at both tails already suffices. This closes the apparent large-C2-error
escape left open by the purely local 4/7 bound **within this parent class**.

### Quantitative three-slice version

Keep the exact target a=1+(t/tau)^2. Suppose at t=(-L*tau,0,L*tau), L>0,
both normalized tensor values are within epsilon<1 of 1. Then

```
1-epsilon <= C <= 1+epsilon,
a_E(+-L*tau)^2 >= (1+L^2)^2*(1-epsilon),
a_E(0)^2 <= 1+epsilon.
```

Both endpoints exceed the middle if
`(1+L^2)^2*(1-epsilon)>1+epsilon`. The concavity lemma therefore gives the
necessary lower bound

```
epsilon >= [(1+L^2)^2-1]/[(1+L^2)^2+1].
```

At L=1 this is 3/5; at L=2 it is 12/13. Equality is inconclusive and is
not asserted to be attainable by a parent. As L grows the bound tends to
1, ruling out any uniform tensor error strictly less than 100 percent for
the all-time target in this ansatz. Only three tensor-value bounds enter;
no small derivatives or symmetric Einstein proper-time intervals are assumed.
Geometry errors require a new scale-factor budget; this is not such a bound.

The explicit q=exp(-3u^2) control has a_E=(1+u^2)exp(-3u^2/2) tending to
zero on both ends, not infinity. Its Einstein Hubble rate decreases, so it
does not contradict the theorem. It demonstrates why the tensor-tail premise
matters. These arguments do not impose NEC on a quantum expectation stress
or on an arbitrary higher-derivative scalar action.

## 3. Exact map to quartic Horndeski: derive its required Jacobian

Could a non-Einstein quartic Horndeski seed instead generate the exact tube?
Use the (+---), X=(partial phi)^2 convention. Write C=C(phi,X), B=B(phi,X),
W=C+BX and

```
X_E=X/W,  dX_E/dX=J_map/W^2,
J_map=W-X*W_X=C-X*C_X-X^2*B_X.
```

For an invertible map J_map must not vanish. In the monotone clock ADM
chart, let V=n^mu partial_mu sqrt(X). The seed extrinsic curvature is

```
K_E^i_j=(K^i_j+s*delta^i_j)/sqrt(W),
s=sqrt(X)*C_X*V/C+sqrt(X)*C_phi/(2C).
```

Quartic Horndeski has no lapse-velocity mixing after its covariant boundary
is removed. Its pulled-back K^2-Kij^2 coefficient is fixed to -M^2/2 by
D's G_T. Expanding `(K+3s)^2-(Kij^2+2sK+3s^2)` therefore gives a K*V
coefficient `-2M^2*sqrt(X)*C_X/C`. The independently derived D ADM action
has coefficient `X^(3/2)*A3`. Hence every such regular map must satisfy

```
C_X/C=-X*A3/(2M^2).
```

The seed Horndeski covariant identity is
`G_T,E=F_T,E-2X_E*(F_T,E)_X_E`. It is replayed against the earlier ADM
coefficients with A1=2F2_X and A3=0. From the inverse tensor dictionary,
D's constant tensor functions throughout the tube require

```
F_T,E=M^2/sqrt(CW),  G_T,E=M^2*sqrt(W)/C^(3/2).
```

Differentiate at fixed phi, using dX_E/dX=J_map/W^2. Multiplying the seed
identity by the nonzero denominators reduces it to
`W-C-X*W_X-X*C_X=0`, or

```
B_X=-2*C_X/X,
J_map=C+X*C_X=C*[1-X^2*A3/(2M^2)].
```

This is a necessary condition derived from the action and tensor map, not
a transcription of a literature transformation formula. Lower-derivative
F/K terms cannot supply K*V or alter the tensor principal identity. The
derivation applies on regular map domains; it cannot license division by
J_map at a zero.

For the pinned D function `A3=4M^2/(1+u^2)^3`, the clock X=1 gives

```
J_map/C=1-2/(1+u^2)^3=Lambda_D/M^2.
u_minus=-sqrt(2^(1/3)-1),  u_plus=+sqrt(2^(1/3)-1).
```

These are exactly two finite zeros: `(1+u^2)^3-2` is even, is -1 at zero,
is 6 at u=1, and has strictly positive derivative 6u(1+u^2)^2 for u>0.
For every finite, positive C, the map loses invertibility at both times.
An infinite/singular C is not an allowed repair. Smooth invertible clock
relabeling cannot remove a zero of the field-map Jacobian. Thus **there is
no everywhere-regular map of the stipulated kind from this exact open D
tube to quartic Horndeski**. This excludes neither approximate deformations
of the tube nor a UV parent with additional fields/operators.

### Two controls against a coordinate overclaim

An illustrative solution of the necessary differential conditions is

```
C=exp[(1-X^2)/(1+u^2)^3],
B=4/(1+u^2)^3 * integral_1^X exp[(1-y^2)/(1+u^2)^3] dy.
```

Here C=W=1 at X=1, even at the two J_map zeros. The metric determinants
stay nonzero on that trajectory, although the change of fields is not
invertible. This is a solution of necessary map equations, not a claim of
a globally valid reconstructed Horndeski theory.

S5's actual auxiliary spatial map has precisely this C but keeps W=1.
Then J_map=1: it is regular, as already proved. It does **not** impose the
additional Horndeski tensor identity; at u=0,X=1 its residual
`1-C-X*C_X` is 2. Removing lapse velocities is not the same as obtaining
Horndeski or ordinary Einstein scalars. The present obstruction leaves
S5's physical reduction and all original crossing certificates unchanged.

## 4. Research consequence

A regular frame change does not rescue the tested conventional parent, and
mapping the exact D tube into quartic Horndeski fails globally. A genuinely
different matching model must introduce new leading dynamics, not merely
rename fields in one of these seed classes. Additional derivative/curvature
operators or fields, tube-deforming corrections with quantitative errors,
and quantum mechanisms remain research questions. No numerical Regge bound,
UV completion, M1 interaction control or full P8 completion is earned.
