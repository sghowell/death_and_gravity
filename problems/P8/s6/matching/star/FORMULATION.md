# P8 S6.16.STAR — constant pairwise star interactions cannot produce a regular flat bounce

This is a bounded background exclusion for a named enlargement of the S6.13
auxiliary matching action. It is not a UV-positivity theorem, a perturbation
health theorem, or closure of the S6 contract or P8(b).

## 1. Action and fixed physical frame

In four dimensions let the physical metric be `h`, with finitely many leaf
metrics `g_i`. The action is

```
S = -G_u/2 int sqrt|h| R_B[h]
    -sum_i G_i/2 int sqrt|g_i| R_B[g_i]
    -2 sum_i int sqrt|g_i| sum_(n=0)^4 beta_in e_n(sqrt(g_i^-1 h))
    +S_m[h,psi].
```

`G_u>=0`, all leaf `G_i>0`, and every Einstein coefficient and `beta_in` is a
finite constant. The beta coefficients have arbitrary signs. Endpoint
coefficients beta0 and beta4 include independent cosmological terms; a
separate central endpoint is equivalent to changing their sum. The source
couples only to `h` and is homogeneous/isotropic on the background, obeys its
equations, and has physical null density `n_h=rho_m+p_m>=0`. Positive-metric
canonical scalar matter is sufficient, not necessary. A bare factor epsilon
would be included in `n_h=epsilon(rho+p)` before canonical normalization; it
cannot suppress a fixed canonically normalized stress by a variable change.

This uses P8(b)'s `+---`, `R_B=-6(DH+2H^2)` convention. The Einstein lapse
equation is `3G H^2=rho`; the null equation is `-2G DH=rho+p`. This is not the
opposite P8(a)/FK curvature convention. The written proof derives the
interaction signs from the displayed action rather than importing them from
a source's beta normalization.

The metrics have a common smooth regular spatially flat FLRW chart on a
connected open interval. All lapse and scale factors are strictly positive
and finite there. Use the positive square-root branch
`sqrt(g_i^-1 h)=diag(N_i,R_i,R_i,R_i)`, where

```
dT=n_u dt,  d tau_i=n_i dt,
R_i=a_u/a_i>0,  N_i=n_u/n_i>0,  c_i=R_i/N_i>0,
H_u=D_T log(a_u),  H_i=D_tau_i log(a_i).
```

No condition `c_i<=1` is imposed. Common-root regularity is an assumption,
not a conclusion about a singular or sign-changing continuation.

## 2. The result and its necessary exception

Call a link **genuine** if `(beta_i1,beta_i2,beta_i3)` is not the zero tuple.
Assume either `G_u>0` or there is at least one genuine link. Every smooth
solution of the full action in the above domain then satisfies

```
H_u(T0)<=0  implies  H_u(T)<=0 for every later T in the same interval.
```

Thus no physical contraction-to-expansion transition is possible, including
a degenerate transition with `H_u'=0` at its first zero. The proof covers
arbitrary zero sets of each `J_i(R_i)=beta_i1+2 beta_i2 R_i+beta_i3 R_i^2`.
It requires neither positive links, a flat vacuum, subluminal cones, an
auxiliary TT inverse nor a constant branch assignment.

If `G_u=0` and every link is endpoint-only, the assertion is not made. With
all beta zero, zero matter and a Minkowski leaf, the physical scale factor
`a_u=1+T^2` is undetermined and bounces. That exception has no central
gravitational equation or viable physical-graviton conclusion. If instead
`G_u>0`, the same disconnected model obeys the usual flat GR null equation.

## 3. Dynamic and algebraic branches are not interchangeable

The exact Bianchi identity is `6 N_i J_i(R_i)(R_i H_u-H_i)=0`.
On a fixed dynamic stratum `H_i=R_i H_u`, with algebraic legs omitted, set

```
K_D=G_u+sum_(i dynamic) G_i/R_i^2.
K_D H_u' - H_u K_D'/2 = -n_h/2.
(H_u/sqrt(K_D))' = -n_h/(2 K_D^(3/2))   if K_D>0.
```

A persistent genuine root has constant `R_i`, constant `H_i`, and
`H_u=c_i H_i/R_i`; its sign is fixed by the positive clock. Including that
leaf in `K_D` is generally false. The actual canonical-scalar solution in
`branches.algebraic_control()` is an explicit counterexample to that wrong
identity, not a counterexample to the theorem.

For the global statement, the proof instead establishes a compact-interval
inequality for the positive part of `H_u`, using bounded logarithmic rates
of the positive ratios and cones. No globally exact all-leg `K` is claimed.

## 4. Evidence and exclusions

The certificate replays literal symbolic lapse/scale and independent exact
Fraction/first-jet calculations, actual algebraic/disconnected controls,
strict input-domain and omission controls, a separately authored covariant
audit, and pinned S6.14/S6.13/adopted-S6 inputs. The branch-switching and
Gronwall argument is a written proof, not a proof-assistant formalization or
a sampled numerical demonstration of every possible zero set.

Excluded: spatial curvature; leaf or multi-metric matter; cycles or general
tensor graphs; variable/scalar-dependent beta or Einstein coefficients;
derivative mixing; negative kinetic coefficients; nonpositive or singular
root/lapse branches; physical NEC violation; and higher-operator/quantum
corrections to this exact action. Their exclusion here is not a no-go for
those models. No full perturbative ghost claim follows from using the HR
pairwise potential. No numerical cutoff, Regge/loop error, exact matching to
the original C/D rows, or universal UV verdict is supplied.
