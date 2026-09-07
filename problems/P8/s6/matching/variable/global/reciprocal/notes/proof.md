# Reciprocal geometry and finite second-metric affine length

This is the independently audited written proof of P8-S6.25.RECIPROCAL.
All hypotheses are in the sibling formulation. In particular
the conclusion is about the specified two-metric domain, not every possible
extension or every notion of completeness relevant to original P8.

## 1. Actual source equation and the center sign

The literal interaction density with all five variable beta coefficients
is the one pinned by S6.20. Its independent f-lapse and f-scale variations
give `rho_If+p_If=(c-y)P/(c*y³)`, where
`P=2(beta1+2 beta2*y+beta3*y²)`. These variations hold the clocks fixed;
there is no erroneous assumption that beta is constant along the solution.
No f matter null stress is added in this theorem.
The beta coefficients in these dimensionless equations mean
`tau² beta_physical/M²`, so the common positive Einstein coefficient is
one. With an unnormalized f Einstein coefficient F, (1) instead has 2F
in its denominator; the sign and affine-length arguments are unchanged.

Since b=alpha/a, `b'=-h*b`. The f Hubble rate and proper-time derivative
are therefore `H_f=b'/(c*b)=-h/c=-z` and `D_f H_f=-z'/c`.
The f metric equation gives the exact undivided identity

```
z' = (c-y)P/(2y³).                                      (1)
```

At the bounce z(0)=0 and z'(0)=h'(0)/c(0)>0. Because P and y are positive,
(1) forces c(0)>y(0). This is not a division by h at the bounce. It also
shows that P=0 or c=y at the bounce is incompatible with the stipulated
nondegenerate smooth reciprocal bounce and source equation.

## 2. The lapse-order inequality cannot reverse

Let `Z=h/y=a*a'/alpha`. Then `Z(0)=0`, and
`B=Z'=(a*a')'/alpha>=0`, with B(0)>0. Hence Z>0 and h>0 for every u>0.
The positive lapse gives z=h/c>0 on this open half-line. Define D=Z-z.
The center sign implies D>0 on some sufficiently short punctured interval.

For u>0, `c-y=yD/z`. Using (1),

```
D' + F D = B,            F=P/(2y²z)>0.                 (2)
```

For any u0 in that initial interval, F is continuous on every finite
interval [u0,u]. The integrating-factor identity is

```
D(u)=exp(-integral_u0^u F)
     [D(u0)+integral_u0^u exp(integral_u0^s F) B(s) ds] > 0.
```

Thus D never vanishes at finite positive time, even when B has zeros.
It follows that c>y and z'>0 everywhere on the future half-line.
There is no assumption that z has a finite limit and no replacement of
the actual lapse by its asymptotic value.

The same geometric hypothesis forces a to diverge. For any sufficiently
small positive u0, `a(u0)*a'(u0)>0`; monotonicity of a*a' gives
`(a²)'=2a*a'>=2a(u0)*a'(u0)>0` afterward. Therefore a(u)->infinity.
It also follows that a>=a(0) on the whole half-line.

## 3. An explicit affine-length bound

For any fixed u0>0, z(u)>=z(u0)>0 afterward, so

```
b*c=alpha*h/(a*z) <= alpha*h/[a*z(u0)].
```

The right side is exactly `-d_u{alpha/[z(u0)*a(u)]}`. Integrating to v and
then taking v->infinity proves

```
integral_u0^infinity b*c du <= alpha/[z(u0)*a(u0)].       (3)
```

No unproved tail expansion is used. A finite interval between the center
and u0 contributes a finite amount because both metric coefficients are
continuous and strictly positive there.

For clarity, the geodesic statement follows directly from the homogeneous
metric rather than a borrowed singularity theorem. For an f null geodesic,
spatial translation invariance conserves the nonzero momentum vector
`p=b² dx/dlambda`. The null constraint gives
`du/dlambda=|p|/(b*c)`, hence `dlambda/du=b*c/|p|`.
Equation (3) bounds its entire remaining affine length. A unit timelike
geodesic with nonzero conserved p obeys

```
ds/du = c/sqrt(1+|p|²/b²) <= b*c/|p|.
```

It too has finite remaining proper length. No conclusion about the
comoving f timelike proper integral is needed or inferred.

For g, the corresponding integrands are a/|p| for null geodesics and
`1/sqrt(1+|p|²/a²)` for timelike geodesics. Since a>=a(0)>0 and u ranges to
infinity, both integrals diverge. Together with spatial homogeneity on
R³ this proves future causal completeness of g. Even/time-reflected
hypotheses prove the analogous past statements.

## 4. Scope and application to the CD-shaped geometry

For `a=(1+u²)²`, alpha=2, the exact formula is
`B=2(1+u²)²(1+7u²)>0`. Consequently every regular lapse and variable-beta
solution in the stated class has the same f affine-length obstruction;
choosing a different lapse within this class cannot remove it. G1's
specific radical reconstruction is not an input to the proof.

The positive relative-shift coefficient is proportional to P, with a
strictly positive prefactor in the actual metric chart. This is why the
positive-P branch is relevant to the previously studied parent. With all
beta coefficients present the TT spring need not equal yP; the theorem
does not assume that beta1-only identity or any scalar health verdict.
The two Einstein tensor cones have relative speed c/y>1 here, but their
high-frequency comparison does not compute an EFT cutoff or establish a
finite-band light-only exclusion.

The result is not an added requirement that original P8 must make both
metrics complete. Nor is a finite affine interval alone a curvature
singularity or a proof that a simultaneous extension is impossible. It
narrows a concrete family of global-parent continuations while preserving
those distinctions. Nonreciprocal scales, f-sector stress, sign-changing
P and other geometric profiles lie outside the theorem, not in a
certified excluded set.
