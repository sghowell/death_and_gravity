# Reciprocal-scale geometric obstruction: specified hypothesis class

P8-S6.25.RECIPROCAL is an isolated conditional theorem, not a new P8 closure
criterion. It generalizes one limitation of exploratory G1 without importing
that candidate's scalar health or modifying any frozen ancestor.

Use the same positive equal Einstein coefficients and physical g frame as
S6.20, in dimensionless g proper time u. On the specified future domain
`[0,infinity) x R^3`, assume the conditions below. Every beta in the
displayed formulas is normalized as `tau² beta_physical/M²`; the positive
common Einstein coefficient has thus been divided out explicitly.
The displayed metrics are `du²-a²dX²` and `c²du²-b²dX²` after rescaling
coordinates and proper lengths by tau. Restoring positive constant units
does not affect the completeness statements.

1. `a` is positive and C2, `h=a'/a`, `h(0)=0`, `h'(0)>0`, and
   `(a*a')'>=0` for u>=0.
2. `b=alpha/a`, with constant alpha>0; `y=b/a=alpha/a²`.
3. `c=N_f` is positive C1 and finite at every finite u, including zero.
4. There is no additional f-sector null stress. The actual f metric equation
   is `-2 D_f H_f=(c-y)P/(c*y³)` in the retained normalization.
5. `P=2(beta1+2 beta2*y+beta3*y²)` is finite, continuous and strictly
   positive. The beta coefficients may depend on the dynamical g clock.

The conclusion is `c>y`, `z=-H_f=h/c` strictly increasing from
zero, and, for every u0>0,

```
integral_u0^infinity b*c du <= alpha/[z(u0)*a(u0)] < infinity.
```

Thus f null and noncomoving timelike geodesics are future incomplete on
this specified domain. The physical g metric is future causally complete.
Time-reflected hypotheses give the corresponding past conclusions. The
quartic CD scale factor `a=(1+u²)²`, alpha=2 obeys the geometric hypothesis.

The theorem does not impose second-metric completeness on original P8,
construct or rule out a simultaneous extension, establish a curvature
singularity, exclude a physical-g bounce, or decide a cutoff/UV question.
It makes no claim after changing the reciprocal scales, adding f null
stress, allowing P to vanish/change sign, or dropping the geometric sign
condition. P positivity is the literal relative-shift sign assumption;
it is not a claim of full perturbative health.
