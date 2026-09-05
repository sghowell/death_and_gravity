# Uniform finite-time hard-channel tree bound

## 1. Scope and local chart

The result is exactly [the operational contract](../FORMULATION.md), not
inclusive scattering or an all-orders theorem. It is uniform in the centre
time, including both infinite tails, on a local hard band and fixed-total-
momentum Hilbert-space fibers on R^3.
All scalar/tensor external assignments and internal polarizations are
included. This is exact arithmetic plus the written lemmas below, with
adversarial self-review, not external review or a kernel-checked proof.

The frozen symbols in `amplitude.py` are separate regression calculations.
This proof does not use them or a frozen propagator: it bounds the exact
time-dependent interaction Hamiltonians and free modes.

Put `delta=1/100`, `ell0=ell(t0)`. The background identities are

```
ell=tau sqrt(1+u^2), ell_dot=x, x_dot=(1-x^2)/ell, H=2x/ell.
```

For `|t-t0|<=delta ell0`,

```
1-delta<=ell/ell0<=1+delta,
|x-x0|<=delta/(1-delta)=1/99,
|log(a/a0)|<=2delta/(1-delta)=2/99.
```

Using `exp(v)<=1/(1-v)` for `0<=v<1`, one gets `1/2<a/a0<2`.
Choose gamma for `|x0|<=9/50`, unitary otherwise. Exact rational checks
give `(9/50+1/99)^2<1/17` and `(9/50-1/99)^2>1/65`. The entire window
therefore stays in the selected previously certified chart. It is not
necessary to stitch these choices into one global vacuum.

Set `L=10^8`, `U=10^9`. External initial |k| is in [L,U], proper subsets
have norm >=L, and a quartic internal wavevector has norm <=2U. Over the
window physical momenta in fixed ell0 units are in [L/2,4U]. Thus

```
q=ell^2*k_physical^2 >=L^2/16>10^14,
q<(8U)^2<Qmax=10^20.
```

Every degree-at-most-four Fourier-mask derivative is bounded by
`Kmax=10^10`. Every required York inverse has nonzero k with k^2>1.
These bounds cover all allowed angles, not just the report's example waves.

## 2. Exact free oscillator estimate

S5.3a supplies `q/2<=Omega^2 ell^2<=3q/2` and
`|partial_t Omega^2|/Omega^2<=10/ell`. For an exact complex free mode u,

```
E=|udot|^2+Omega^2|u|^2,
Edot=(partial_t Omega^2)|u|^2,   |Edot|<=10E/ell.
```

The full window has length `2delta ell0`, so
`E(t)/E(left)<=exp(20/99)<=99/79<2`. In fixed ell0 units the prescribed
initial data give `E(left)=Omega(left)`. Let `W=10^11`. The band implies
`Omega(left)<5U<W/4` and `Omega^2>L^2/16`. Hence, throughout I,

```
E<W/2,   |u|^2<8W/L^2<1,   |udot|^2<W/2<W^2.
```

Each normalized oscillator leg, including each end of an internal
propagator, therefore has `|Y|<=1`, `|P|<=W`, also for conjugate modes.
This energy estimate is exact, not an error bound for a WKB approximation.
The curvature scale is at most `sqrt(6)/ell`, far below the chosen band.

## 3. Canonical seeds and fixed-unit coefficients

Write `s=ell/ell0`, normalize a0=1. Compact scalar K obeys
`1/8<K<=65/2` on the covering domains. Fixed-unit K is unchanged in unitary
and divided by s^2 in gamma. A common weaker bound is `1/32<Kfixed<130`.
The exact linear canonical map is

```
Q=Y/[a^(3/2)sqrt(2Kfixed)],
p_chart=[2Kfixed(P-fY)+beta Y]/[a^(3/2)sqrt(2Kfixed)],
beta=-2Kfixed Bphase,   f=3H/2+dot Kfixed/(2Kfixed).
```

p_chart is momentum density divided by a^3. This is an exact canonical
map with a quadratic generating boundary, not substitution of the full
interacting equation P=Ydot. P=Ydot is valid in the free interaction
picture. All additional generating terms affect H2, not H3/H4.

The compact f majorants are 756 (unitary) and
`2460299464273/21474836480` (gamma), both <`Fmax=1000`. Fixed-unit f is
bounded by 2Fmax. Further,

```
a^(-3/2)/sqrt(2Kfixed)<12,
|Bphase|<=8Qmax+4,    2Kfixed<=260.
```

Thus |p_chart| is at most `12[260(W+2Fmax)+260(8Qmax+4)]`. The gamma
canonical swap is `zeta=p_gamma/(2q_fixed)`, `p_zeta=-2q_fixed b`;
q_fixed>=1 bounds zeta by the previous momentum majorant and
`|p_zeta|<=96Qmax` suffices. Each entry of a normalized TT polarization
has absolute value <=1. Tensor coordinate and momentum entries are bounded
by 6 and `2(W+6)`. Every seed is therefore <the common **B=10^30**.
The allowed a,s rescalings and the normalization derivative are accounted
for explicitly.

For the invariant monomial `sigma^p rho^r S2^w`, the fixed-unit coefficient
is the compact S5.1 coefficient times `s^(p+2r+2w-2)`. At weighted order
<=4 this exponent is in [-2,6], so its absolute unit-conversion factor is
<=64 on `s in [1/2,2]`. The compact bounds come from the pinned S5.1
certificate; no unproved momentum-uniform continuity claim replaces them.

## 4. Positive-series majorant induction

`Series` stores nonnegative norms through degree four. The degree-n norm
is the sum of absolute values of all labelled-mask coefficients of that
degree. Products are bounded by nonnegative convolution, derivatives by
Kmax multiplication. Four possible labels per field also bound cubic cases.

Let `b=4B e`, where e marks field order, and `r_g=3b` bound each entry of
`2zeta I+gamma`. Coefficientwise bounds are

```
g: 1+r_g,
g_inverse: 1+sum_(n=1)^4 3^(n-1)r_g^n,
det(g)-1: positive-degree part of 6(1+r_g)^3.
```

The factors count matrix products and determinant permutations. The exact
determinant background remains 1, not the overestimated constant 6.
Finite binomial series with absolute coefficients, exponents 1/2 and -1/2,
bound the volume and inverse volume. Full coordinate geometry gives

```
Gamma <= (9/2) Kmax g_inverse r_g,
R <= 9 g_inverse [6Kmax Gamma+18Gamma^2].
```

The free momentum bound is `pi_free<=4+(4+1/6+1+4)b`, since |H|<=4.
At n=1,2,3 take the degree-n source from `3Kmax pi+9Gamma pi`.
The exact York inverse implies `|W|<=2|source|/k_min^2<=2|source|`, and
`|LW|<=4Kmax|W|`; add `8Kmax source_n` at that degree. This includes all
nonlinear vector constraints. The prior trace lemma proves W4 cannot enter
H4, so no uncomputed fourth-order inverse is assumed.

Set `pMixed<=3pi g/sqrt(g)`. Bound sigma by the positive-degree part of
`3pMixed`, and S2 by the degree-at-least-two part of `12pMixed^2`. Their
removed background/linear parts vanish by exact prior identities; arbitrary
terms are not being discarded. Insert these and R in every stationary
monomial with its fixed-unit coefficient bound. Multiply by the volume,
retain `-2H pi:g` with majorant `72pi g`, and multiply by 8 for a^3<=8.
The explicit Adot and gamma generators are at most quadratic and do not
contribute to H3/H4 in this linear metric/canonical parametrization.

The resulting positive series h bounds all allowed physical scalar/tensor
leg assignments. A redundant `6!=720` factor also covers normalized
repeated-mode occupation/Wick factors for four external legs and two ends
of one internal contraction. Set `B3=720h_3`, `B4=720h_4` and record their
exact integers. These bounds follow from finite sums and the generic inverse,
not from sampled vertices. Tests independently check the positive-series
algebra against polynomial/binomial expansion and every scale inequality.

## 5. Exact interaction-picture contacts and exchange

The verified physical scale restoration gives, in fixed patch units,

```
H_int=H3/(M ell0)+H4/(M ell0)^2+higher field orders.
```

Use exact free modes from section 2. At a quartic tree there are three
external 2+2 partitions and three internal polarizations. Bounding both
time orderings separately gives the safe factor 18. The dimensionless
window length is 1/50<1. Each retained matrix entry is therefore at most

```
cubic:           B3/(M ell0),
quartic tree:    [B4+18B3^2]/(M ell0)^2.
```

The second line includes one H4 contact and two H3 insertions, with scalar
and both tensor propagation. In velocity variables the quartic Legendre
contact and differentiated-propagator contact terms must be combined
consistently. Here canonical time ordering acts directly on Y and P; no
spurious delta term is added by differentiating an already time-ordered
propagator. Quantum ordering effects first contributing loops are excluded
by the tree contract, not claimed small.

Use Fourier measure `d^3k/(2pi)^3` and ordered momentum coordinates followed
by normalized Bose symmetrization. The latter's finite combinatorial factors
are covered by 720. Spatial homogeneity gives an exact total-momentum delta.
Decompose the one- and two-particle spaces into total-momentum fibers before
taking an operator norm; the coordinate change `(k1,k2)->(P=k1+k2,k1)` has
Jacobian one. No delta distribution is bounded in absolute value.

At fixed P a one-particle fiber has three species labels. A two-particle
fiber has one relative momentum and at most nine species pairs, so its
measure is at most `9(4U)^3/(2pi)^3<3N`, where `N=3(4U+1)^3` is a loose
momentum/species **measure bound**, not a physical state count. Row and column
integrals of the bounded kernels are therefore at most D times the entrywise
bound, with the deliberately oversized common `D=4N^3>3+3N`. The Schur
test bounds both rectangular cubic 1-to-2/2-to-1 blocks and quartic 2-to-2
blocks by D times their respective entrywise bounds. Uniformity in P then
bounds the corresponding direct-integral operators on normalized wavepackets.

For quartic transitions apply the hard-transfer mask first; the removed
forward entries are not bounded. Momentum conservation fixes the one internal
tree wavevector, so there is no loop integral hidden in 18. Vacuum-to-three
production is deliberately excluded: its total-momentum delta cannot be
treated as a bounded map from the vacuum to an infinite-volume Fock sector.
No claim about a gravitating torus with unresolved homogeneous constraints
is needed. The normalizations are the continuum canonical ones; no box
volume factor is inserted or sent to infinity.

Set `C3=D B3`, `C4=D(B4+18B3^2)`. Since ell0>=tau, it suffices to impose

```
M tau>=1000 C3,    (M tau)^2>=1000 C4.
```

The first power of ten meeting both exact rational inequalities is
**10^270**. One common finite choice works at every t0, above curvature
on the specified hard band. It is not a necessary scale, an actual cutoff
estimate or a realistic bounce-duration proposal. Cancellations of frozen
vertices cannot invalidate a sufficient absolute bound of this kind.

## 6. What remains unproved

Removed soft/forward channels, homogeneous backreaction, vacuum production,
full Fock-space norms, loops, higher tree orders, arbitrary
occupation numbers and accumulated evolution over all time are not bounded.
No global vacuum or chart-independent in/out amplitude is constructed.
This is the scaled D-only M0 family, not M1. The selected tree criterion
is met; unrestricted EFT/UV completion and the broader P8 program remain
separate targets.
