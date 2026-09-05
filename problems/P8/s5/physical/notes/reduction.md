# Exact physical-mode reduction through fourth order

## 1. Local units and input

Use the S5.1 D-only M0 Hamiltonian, with `M=ell=a0=1` at the evaluation
point. Here `ell=tau sqrt(d)`, `d=1+u^2`, `x=u/sqrt(d)` and `z=1-x^2`.
This is a constant change of units in each local patch; it is not a global
running canonical normalization. The compact background functions are

```
P=1-6x^4+10x^6-4x^8,    1/4<P<=1,
J=2P,  Theta=H=2x,  Lambda=1-2(1-x^2)^3,
Hdot=2-4x^2,    Hdot+3H^2=2+8x^2.
```

The stationary S5.1 density `h(sigma,rho,S2)` retains invariant degrees
0--4, with weights `(1,1,2)`. Its lapse equation is solved through degree
three, including the second-order correction needed in the quartic density.
It has no inverse Theta or Lambda. S5.1 supplies its all-time denominator
proof and the background equations; neither is redefined here.

## 2. Spatial gauge and exact canonical coordinates

Choose the linear metric parametrization

```
hat_h_ij=a^2 g_ij,   g_ij=(1+2 zeta)delta_ij+gamma_ij,
partial_j gamma_ij=0,   gamma_ii=0.
```

This is a valid *small-field* spatial gauge. At the background its
gauge-fixing derivative on a spatial diffeomorphism is the flat York
operator of section 3, invertible at every nonzero spatial momentum.
Order-by-order gauge fixing uses this same linear inverse. The claim is
a finite formal perturbative construction, not a global gauge or an
existence theorem on an unrestricted Sobolev space. In particular, large
fields may make `g` nonpositive and lie outside this patch.

Let `pi_bar=a^2 pi`, `A=a^3 H`, and use flat contractions throughout this
section. Decompose the contravariant momentum **density** as

```
pi_bar = [-A(1+zeta)+p/6]I + Pi_TT + A gamma + L W,
(LW)_ij=partial_i W_j+partial_j W_i-(2/3)delta_ij partial_k W_k.
```

The longitudinal piece is flat-tracefree and orthogonal, after spatial
integration by parts, to `dot gamma_TT`. Thus `p` and `Pi_TT` are exact
canonical momenta for this gauge, not merely their linear approximations.
Writing `gamma:gamma=sum_ij gamma_ij^2`, the one-form obeys

```
integral pi_bar:dot g
 = integral [p dot zeta+Pi_TT:dot gamma] + dot B
   + Adot integral [6 zeta+3 zeta^2-gamma:gamma/2],
B=integral [-6A zeta-3A zeta^2+A gamma:gamma/2].
```

The time-dependent relation `hat_h=a^2 g` supplies an additional
`2H pi_bar:g` in the symplectic Lagrangian. Consequently the reduced
Hamiltonian, before the gamma swap, is

```
Hred = sqrt(g) h(sigma,R[g],S2)
       -2H pi_bar:g -Adot(6 zeta+3 zeta^2-gamma:gamma/2).
```

These boundaries are only linear/quadratic in this parametrization but
are essential: dropping them gives incorrect tadpoles or the old quadratic
Hamiltonian. No nonlinear Darboux correction is hidden in the code.

## 3. All three momentum constraints

For a contravariant density of weight one the connection traces in its
divergence cancel. The spatial constraints reduce to

```
C_i=partial_j pi_bar^{ij}+Gamma^i_jk(g) pi_bar^{jk}=0.
```

The sign and overall factor relating this to the ADM shift constraint
do not affect its zero set. The exact flat operator on `W` and inverse are

```
E=Delta I+(1/3)grad div,
E(k)=-k^2 I-k k^T/3,
E^-1(k)=-(I-k k^T/(4k^2))/k^2,
det E=-(4/3)(k^2)^3.
```

At degree n insert the solutions through n-1 into `C`, take its degree-n
source, and set `W_n=-E^-1 source_n`. All additional coefficients multiplying
`W_n` have positive field degree, so this kills precisely the degree-n
residual. Induction proves the finite recursion for arbitrary nonexceptional
momenta. The leading result is `W_1=-grad Delta^-1 p/8`. Both `W_2` and
`W_3` are generally nonzero. Already two scalar inputs generate a transverse
part of `W_2`; a scalar-only shift ansatz would miss interactions.

Through fourth Hamiltonian order only `W_1,W_2,W_3` are required. A possible
`LW_4` contribution multiplies the Hamiltonian's isotropic background first
variation and vanishes by flat tracelessness. All other terms containing it
have degree at least five. The same argument applies to the explicit
`-2H pi_bar:g` term. This is why the full-degree-zero integrated wavevector
never needs to be inverted in a nonexceptional n-leg kernel.

The implementation rejects a zero proper-subset wavevector before
construction. The solver also raises on an actual nonzero homogeneous
constraint source. It does not claim homogeneous backreaction is absent.

## 4. Geometry, invariant substitution and labelled Fourier jets

Compute `g^-1`, `sqrt(det g)`, the full coordinate Christoffels and `R[g]`
by formal finite Taylor arithmetic. Define

```
pMixed=pi_bar g/sqrt(det g),
sigma=tr(pMixed)+3H,
S2=tr(pMixed^2)-tr(pMixed)^2/3,    rho=R[g].
```

Insert these into the stationary density and multiply by the **full** volume
factor before extracting the target order. No conformal-only curvature
shortcut or frozen determinant is used. As independent geometry checks,
a pure coordinate pullback of the flat metric has `R=0` through fourth
order, and `g=(1+2 zeta)I` gives

```
R=-4 Delta zeta/(1+2 zeta)^2+6(grad zeta)^2/(1+2 zeta)^3.
```

For n external legs work in the ring generated by labelled amplitudes
`epsilon_i` with `epsilon_i^2=0`, carrying Fourier momenta `k_i`.
A bitmask selects a product of distinct labels. Multiplication convolves
only disjoint masks, and a derivative multiplies by `i sum_mask k`.
Extraction of the full mask is exactly the mixed derivative with respect
to all n labels. Repeated fields use distinct labels; e.g. the cubic
coefficient of `(epsilon_1+epsilon_2+epsilon_3)^3` is 6, not 1. Finite
binomial/geometric series are exact in this quotient. This establishes
the procedure for any allowed external data without claiming the finite
report fixtures enumerate continuous momentum space.

## 5. Quadratic bridge and the gamma canonical swap

The calculation reproduces, with `q=|k|^2` in fixed local units,

```
H2_scalar=-q zeta^2+(Theta p/2-Lambda q zeta)^2/J,
H2_tensor=2 Pi_TT:Pi_TT+q gamma:gamma/8.
```

There is no tensor momentum-coordinate or scalar-tensor quadratic mixing.
For a polarization `E`, choose `gamma=E Q`, `Pi_TT=E P/(E:E)`; this is a
dual canonical basis, even when the rational polarization has norm other
than one. The two independent TT polarizations span the full tensor space.

Near gamma crossing use the exact linear canonical transformation

```
p=-C b,    Pb=C zeta,    C=2 a^3 k_physical^2.
```

At a fixed comoving k, `Cdot=H C`. Integrating `-C b dot zeta` by parts
therefore adds **`-H b Pb`** to the transformed Hamiltonian. At the local
evaluation point `C=2q`, giving

```
H2_gamma=(q Lambda^2-J)Pb^2/(4Jq)
         +(Theta Lambda q/J-H)b Pb+Theta^2 q^2 b^2/J.
```

The same linear transformation applies to H3 and H4 without further
higher-degree generating terms. The phase Hamiltonians remain regular at
Theta=0 and Lambda=0. A velocity chart additionally needs invertible
`H2_PP`; its finite-frequency coefficient is

```
K_unitary=J/Theta^2,
K_gamma=q J/(q Lambda^2-J).
```

The S5.1 covering intervals give the sufficient domains stated in the
contract. The positive gamma denominator restriction also applies at
internal momenta in section 7. A chart failing that test does not establish
an unhealthy physical mode or a physical cutoff; use phase variables or
another admissible chart.

## 6. Canonical normalization and its physical time derivative

In either velocity chart the linear momentum response is

```
P0=alpha Qdot+beta Q,  alpha=2K,
beta_unitary=2 Lambda q/Theta,
beta_gamma=-alpha(Theta Lambda q/J-H).
```

For tensors `Pi_TT=dot gamma/4`. Unit kinetic normalization is
`Qc=Z Q`, with `Z=sqrt(2K)` for the scalar and `Z=sqrt(E:E)/2` for
a tensor polarization. The scale `M` is restored separately. Crucially,

```
Q=Qc/Z,   Qdot=(dot Qc-(D log Z)Qc)/Z.
```

To evaluate the derivative, a coefficient written `f=d^w f_old` obeys
the fixed-patch physical derivative

```
D_w f=(1-x^2) partial_x f -2x q partial_q f -2w x f,
D log Z=(D_w K)/(2K),
w=0 for unitary,  w=1 for gamma b/ell0.
```

Here `q=k_comoving^2/d` labels compact local momentum. Its derivative in
moving compact units is `-2xq`; the weight term removes the derivative of
the unit itself. This formula follows directly from `ell partial_t x=1-x^2`
and the chain rule at fixed comoving momentum. The gamma formula is checked
against differentiation of the original uncompactified kinetic expression.
The tensor normalization and polarization are constant in the fixed patch.

The code factors the product of external `1/Z` outside the rational jet
calculation; it still inserts `-D log Z` for each normalized coordinate leg.
After this substitution the labelled two-velocity kernel is exactly 1,
corresponding to `dot Qc^2/2`. The background `a^3` measure remains present
in the action; this is not yet a stationary Minkowski scattering frame.

## 7. Quartic Legendre correction and tensor completeness

Write `H=H2+H3+H4+...`, let `A=H2_PP`, and solve the velocity equation as
`P=P0+P1+...`, where `P1=-A^-1 H3_P(P0)`. Stationarity cancels the term
linear in P1 and gives

```
L2=P0 Qdot-H2(P0),
L3=-H3(P0),
L4=-H4(P0)+(1/2)H3_P(P0) A^-1 H3_P(P0).
```

This identity is independently checked with generic scalar polynomial
coefficients. Its matrix/functional version follows by the same symmetric
quadratic completion. There is no derivative of `A` in this pointwise
momentum elimination; its spatial Fourier inverse is evaluated at the
contracted momentum.

For four labelled legs use the three unordered partitions `01|23`,
`02|13`, `03|12`. Compute each source using a cubic Hamiltonian kernel
with the two chosen external legs and one raw momentum leg carrying minus
their sum. Contract two scalar sources with `A_scalar^-1=alpha(K)`.
For each orthogonal TT polarization E at the internal K, contract tensor
sources with `A_E^-1=(E:E)/4`. Both polarizations are summed. This dual
basis rule follows from `H2_tensor=2P_E^2/(E:E)` and is basis independent.

There is no remaining 1/2 when summing the three **unordered** partitions;
each partition's two ordered versions cancel that factor. A toy
`H3=g P Q^2`, `H2=P^2/2` yields the labelled quartic coefficient `12g^2`,
agreeing with `L4=g^2 Q^4/2`. Tensor contributions are nonzero even for
some all-scalar external configurations and are explicitly tested.

This correction is a contact term generated by the Legendre transform.
It is **not** the physical exchange amplitude of two cubic vertices joined
by a propagating scalar or graviton. S5.3 must still compute that amplitude
and combine it with the full L4 contact term before inferring control.

## 8. What the compact-time bounds do prove

For fixed allowed momenta, the York inverses are finite constants. Every
geometry and constraint operation is polynomial in the perturbations,
and the background and canonical boundary coefficients are polynomial in
x. The only compact-time denominators in phase kernels are powers of P
inherited from the stationary lapse solution. Thus the finite-order phase
kernels are rational functions of x with no poles on `[-1,1]`.

The report checks this denominator form for each listed symbolic phase
kernel and bounds it by the numerator's coefficient L1 norm times `4^m`
divided by the constant denominator factor, using `P>1/4`. This proves
boundedness including both tail limits **for those fixed-momentum kernels**.
It does not bound growth as frequency rises or momenta become exceptional,
nor establish cancellations, adiabaticity, a physical amplitude or a cutoff.
