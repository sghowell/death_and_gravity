# G1 feasibility argument — not an independently reviewed certificate

This working note records a new candidate and the limitation found in it.
The exact tests support the algebra below; no certified ledger promotion
or original P8 closure follows. Definitions are in the
[candidate formulation](../FORMULATION.md).

## 1. Reconstruct before choosing the new lapse

Work temporarily in `M=tau=1`, with `a'=ah`, `b=2/a`, `y=2/a^2` and
an arbitrary odd geometric profile z. Put `c=h/z`, so `H_f=-z` and
`D_f H_f=-z'/c`. The literal interaction lapse polynomials are
`U=b0+3b1*y`, `V=b1+b4*y^3`, and `P=2b1`.

For the coefficient reconstruction in the formulation, all four metric
equations are identities:

```
3h^2=n/2+2U,            3z^2=2V/y^3,
-2h'=n+(y-c)P,          2z'/c=(c-y)P/(c*y^3).
```

The exact jet calculation also verifies the separate free-chi equation,
the sourced clock equation multiplied by its positive clock speed,
and the undivided sourced Bianchi relation:

```
w'+3hw=0,              w=1/(10a^3), k=n-w^2,
k'/2+3hk+2[b0'+(c+3y)b1'+c*y^3*b4']=0,
3P(-yz-h)-2(b1'+y^3*b4')=0.
```

The derivative acts on the coefficient profiles, not on y or c held fixed
in a partial scalar derivative of the interaction. In the jet engine,
a,h,h',z,z' are independent inputs and the total derivative retains
`a'=ah`, `h''` and `z''`. No scalar equation is inferred solely from a
metric fit, no unsourced clock equation is imposed, and the optional f
scalar is absent. Isotropy and the diagonal root make the remaining
background vector/off-diagonal equations vanish. No perturbative rank
or stability statement is obtained from that background symmetry.

## 2. Analytic positive lapse and separation from c=y

For the actual CD scale factor,
`A=h'/y^3=(1-u^2)(1+u^2)^10/2` is an even polynomial. Define

```
z'=(A+sqrt(A^2+2))/2=1/(sqrt(A^2+2)-A),  z(0)=0.
```

This is strictly positive, even and real analytic for every finite u.
Its integral is odd, analytic and has the sign of u. Taylor expansion gives
`z=u+u^3+O(u^5)`, hence
`c=h/z=4-8u^2+O(u^4)`. The zero in that quotient is removable, not a
degenerate lapse. In particular, this action differs from the constant-c=4
action already at second order in u.

For u>=0 put `Z=h/y=2u(1+u^2)^3` and
`B=Z'=2(1+u^2)^2(1+7u^2)>0`. The positive root z' of
`r^2-Ar-1/2=0` is smaller than B whenever `B(B-A)>1/2`.
With `x=u^2`, this comparison on x in[0,1] is a degree-14 polynomial.
On each of sixteen equal rational x intervals its fifteen Bernstein
coefficients are strictly positive. All 240 coefficients are recomputed
with a separate Fraction/binomial implementation and compared to the
symbolic coefficients; their overall minimum is exactly `5/2`.
The convex-hull property of the Bernstein basis supplies a continuous
interval bound, not a sampling inference.

For u>=1, A<=0 and `z'<=sqrt(2)/2<1<B`, so the same comparison holds
without polynomial subdivision. Integrating from the common zero gives
`0<z<Z` for u>0. Therefore `c=h/z>h/Z=y>0`. Evenness gives the same
strict lapse/root separation for negative u. At zero it holds because
`c=4>2=y`. Every denominator in the reconstructed coefficients is thus
nonzero at every finite time.

## 3. Positive clock on the entire real line

The derivative of A is exactly
`u(1+u^2)^9(9-11u^2)`. Its global maximum is
`(20/11)^10/11<40`; after u=1 it is nonpositive and decreases to minus
infinity. Since `sqrt(A^2+2)+A` is increasing in A,

```
sqrt(A^2+2)-A = 2/(sqrt(A^2+2)+A) > 2/81.
```

Here `sqrt(40^2+2)<41`. Also `y^3/w^2=800` exactly. Consequently
`n/w^2>1600/81` and `k=n-w^2>(1519/81)w^2>0` at every finite time.
This is a bare canonical-clock bound, not the coupled scalar no-ghost
condition. Dropping the smoothing term 2 is a negative control: at the
center the resulting k is `-1/100`.

The positive analytic square root of k gives a locally invertible
analytic scalar clock. Section 4 gives `k~8/u^2` in each tail, so its
integral diverges logarithmically in both directions. Thus phi maps R
onto R monotonically and has a global real-analytic inverse. The stated
beta functions are genuine analytic functions of this dynamical phi
for every finite phi; they are not external functions of coordinate time.

## 4. Completeness, curvature and the lost massive threshold

For sufficiently large u, `-A>=u^22/4` and
`0<z'<=1/(2|A|)<=2/u^22`. The finite positive limit
`z_infinity=integral_0^infinity z'(v)dv` therefore exists, and
`z_infinity-z(u)<=2/(21u^21)`. Rationalizing the derivative gives the
exact leading limits

```
z'~u^-22,       c~4/(z_infinity*u),
n~k~8/u^2,      phi~sqrt(8)*M*log u  (future tail).
```

Odd/even symmetry determines the past. In physical units,
`H_g=h/tau` tends to zero in both tails;
`H_f=-z/tau` tends to the two opposite finite constants. Its proper-clock
derivative `-z'/(c*tau^2)` tends to zero and is analytic through the
center. Both FLRW curvature tensors have bounded orthonormal components
on this all-time physical-g domain. No comparison to a computed cutoff
or quantum correction is being made.

The physical metric g has `a>=1`, so its null affine integrals `integral a dT`
and all timelike proper integrals `integral dT/sqrt(1+p^2/a^2)` diverge in
both directions. **This completeness assertion is only about g.** Although
the f comoving proper time also diverges logarithmically, its null affine
integral `integral b*c du` is finite in either tail. Noncomoving f timelike
geodesics likewise have finite proper length. Bounded curvature does not
turn those f-domain statements into geodesic completeness, nor does this
note construct a simultaneous extension of the full two-metric solution.

The exact relative shift coefficient is positive because `P=2b1>0`.
The two full tensor principal speeds relative to physical g are 1 and
`c/y>1`. Thus this is not an all-parent-modes subluminal construction.
A light-only interpretation would need its own controlled reduction.

The literal TT stiffness is `mu=yP`, giving the dimensionless algebraic
relative mass coefficient

```
tau^2*m_alg^2 = 2*y*z'*(y^3+c)/(c-y) ~ 4/u^30.
```

The tail calculation substitutes the established limit z_infinity only
in the continuous leading coefficient; it does not replace the actual
background by a constant-z solution on a finite interval. The exact
rationalized limits are checked in the module. This coefficient is
positive at finite time but has infimum zero over the complete history.
It cannot serve as a uniform positive massive threshold.

That last fact is not itself a theorem about a rolling spectral gap, a
physical instability, or every possible decoupling mechanism. The kinetic
weights and physical source overlaps also change in the tails. Proving
that these permit or obstruct an effective light theory requires the
full canonical response, state and error estimates. The construction
therefore removes one background-domain obstacle but does not close the
matching or UV questions, or supply an original C/D witness.

## 5. Verification status

The nineteen ordinary exploratory tests pass, including all seven metric/
scalar identities, all 240 Fraction/Bernstein coefficients, the analytic
center, global-clock inequalities and leading rationalized tail limits.
Ruff passes. An independently authored scientific audit, frozen formulation
review and source-hashed report have not been performed. No CERTIFIED
ledger row or claim of original-P8 closure is made for G1.
