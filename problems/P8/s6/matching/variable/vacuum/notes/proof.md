# Local proportional-vacuum obstruction

This is the independently reviewed proof for S6.28.LOCAL_VACUUM. It is a continuous
argument about the frozen S6.20 action on its specified field interval,
not a search on a numerical grid.

## 1. Action, field range, and independent vacuum variables

Use the literal action and positive Einstein coefficients in the
formulation. On the rolling solution S6.20 defines

```
d=1+u², a=d², y=2d^-4, h=4u/d, h'=4(1-u²)/d²,
nbar=2(y³/c-1)h',  kbar=nbar-1/(100a^6),
b1=y³ h'/[c(c-y)], b2=b3=0,
b4=3h²/(2c²)-b1/y³,
b0=3h²/2-nbar/4-3b1y.
```

Here c is a fixed parameter labeling an action, either 1 or in (2,4].
It is not a variable to solve for in a vacuum of that same action.
The strict positive kbar bounds of S6.20 make

```
phi=M integral_0^u sqrt(kbar(v,c)) dv
```

an analytic one-to-one label of the specified closed field interval. Its
inverse defines `beta_n(phi)=M² b_n(u(phi/M),c)/tau²`. We use this label
u when evaluating coefficients at a proposed constant phi*. In that
evaluation u is not a spacetime coordinate or an imposed rolling clock.
In particular, nbar and kbar in this proof are coefficient-profile
identities, not nonzero matter velocities imposed on the vacuum.

Write the stationary coframes as `f=r²g` with r>0. Rescale coordinates to
g=eta if desired. The ratio r is independent of the action label c and
the coefficient-profile value y(u). Both canonical scalars are constant;
the actual canonical kinetic stresses vanish. The free chi equation is
then automatic, with arbitrary constant chi value.

## 2. Full metric and clock algebra

Before making the proportional restriction, the literal interaction
density in four homogeneous coframes is

```
-2[Ng ag³ beta0 + beta1(Nf ag³+3Ng ag² af)
   +3beta2(Nf ag² af+Ng ag af²)
   +beta3(3Nf ag af²+Ng af³) +beta4 Nf af³].
```

These coefficients follow from
`(1+t Nf/Ng)(1+t af/ag)^3`. Differentiate both lapses and both scales,
holding the clock fixed, before setting `af=r ag,Nf=r Ng`. The interaction
densities are

```
rho_g=2U_beta,
rho_f=2V_beta/r³,
U_beta=beta0+3r beta1+3r² beta2+r³ beta3,
V_beta=beta1+3r beta2+3r² beta3+r³ beta4.
```

Both pressures are minus their respective densities. The flat Einstein
tensors vanish. Thus the two lapse conditions `U_beta=V_beta=0` also
enforce the spatial equations; no lapse or source equation is omitted.
The literal canonical clock Euler equation at a constant field is

```
2[beta0,phi+4r beta1,phi+6r² beta2,phi
  +4r³ beta3,phi+r⁴ beta4,phi]=0.
```

The binomial factors are those of `e_n(r I)=binomial(4,n)r^n`, not the
three-dimensional coefficients in either metric lapse equation. The
ratio r is held fixed when varying the clock. Before stationarity its
left side is added to `Box_g phi`; the sign follows from the positive
canonical action and the literal -2 interaction normalization.

For the actual b2=b3=0 profiles the complete algebraic conditions are

```
U_vac=b0+3r b1=0,
V_vac=b1+r³ b4=0,
[2M/(tau² sqrt(kbar))] [b0,u+4r b1,u+r⁴ b4,u]=0.
```

Indeed `dphi/du=M sqrt(kbar)` and
`beta_n,phi=M b_n,u/(tau² sqrt(kbar))`. The prefactor is finite and
strictly positive for every finite admitted action parameter. No
vanishing clock force is inferred from a separate-conservation assumption.

The metric equations alone will already be inconsistent below. The clock
equation is retained nevertheless: their consistency would not by itself
give a vacuum. Generic coefficient-point controls test this distinction.

## 3. The unique f-flat root

The exact interval bounds inherited from S6.20, supplemented by the
immediate h² bound, are

```
0<=h²<=4/25, h'>3, 1<y<=2, y³>7.
```

For example `d<=101/100` gives `d^12<8/7`; all displayed clearing
denominators are positive. Define

```
theta=3h²(c-y)/(2c h'),
b4=-(b1/y³)(1-theta).
```

For c=1, c-y<0, so b1<0, theta<=0 and b4>0. For 2<c<=4,
c-y>0, b1>0, and

```
0<=theta<=3h²/(2h')<2/25<1,
```

so b4<0. Both b1 and b4 are nonzero throughout the field interval.
Consequently the exact f-flat equation has precisely one positive root:

```
r_*³=-b1/b4=y³/(1-theta).
```

For c=1, `0<r_*<=y`. For the positive branch,

```
y<=r_*,  r_*³<200/23<(9/4)³,
```

where the last difference is the strictly positive rational number
`3967/1472`. Thus `r_*<9/4` uniformly. This is a continuous bound,
not a sampled solution of a cubic. The singular action parameter c=2
is never inserted into the profiles.

## 4. A strictly nonzero g residual at that root

For any positive r let `S=r²+ry+y²`. Direct algebra with the actual
profiles gives the exact residual identity

```
U_vac+3y³ V_vac/S
 =3h²/2-nbar/4+9h²y³r³/(2c²S).
```

Impose the f equation exactly, so r=r_* and V_vac=0. Because
`S-3ry=(r-y)²>=0`,

```
U_vac <= (3h²/2)[1+y²r_*²/c²]-nbar/4.
```

For c=1, `r_*<=y<=2` and `nbar>36` imply

```
U_vac < (6/25)(1+16)-9 = -123/25.
```

For 2<c<=4, `r_*<9/4`, `y<=2`, and
`nbar/4>(3/2)(7/c-1)` imply

```
c² U_vac < Q(c),
Q(c)=(87c²-525c+243)/50.
```

The complete continuous lapse interval is controlled by the exact identity

```
-459/50-Q(c)=3(c-2)(117-29c)/50>=0,
```

since `117-29c>=1`. Equivalently, on c=2+2z the margin has the
nonnegative degree-two Bernstein coefficients `[0,177/50,3/25]`.
The zero at c=2 belongs only to this polynomial bounding interval, not
to an admitted action. Dividing by c²<=16 gives

```
U_vac < -459/(50c²) <= -459/800.
```

The actual flat g density would be `2M² U_vac/tau²`, hence is strictly
negative, while the flat g Einstein tensor and constant canonical-matter
kinetic density are zero. The g metric equation cannot hold. This
excludes all simultaneous solutions of the two metric and clock equations
in the stated domain, regardless of whether the clock equation vanishes.

The density gaps have units M²/tau² and apply after imposing the f-flat
equation exactly. They are not a claim that arbitrary independently
bounded approximate metric/source errors obey the same gap. A quantitative
multi-equation error budget would require additional control of the f
residual and the ratio.

## 5. Omission and domain controls

At u=0 the coefficients are even, so every first clock derivative vanishes.
For every admitted c, r_*=2 satisfies the f and clock equations but
`U_vac=2-16/c<0`. At c=4 this is -2. Omitting the g equation therefore
produces a false vacuum. Conversely, c=4, u=0, r=13/6 solves the g and
clock equations but gives `V_vac=-469/432`, exposing the omission of f.

A separately named constant coefficient point with r=2,
`(b0,b1,b2,b3,b4)=(-6,1,0,0,-1/8)` and zero clock derivatives satisfies
all vacuum equations. Keeping those metric values but choosing
`b0,phi=1` makes the clock residual 2. These are controls of the equation
system, not assignments to the frozen coefficient profiles.

Adding a constant 2 to the actual c=4 b0 would also produce a flat
stationary point at phi=0,r=2. That changes the action on the specified
bounce field interval; it is not an allowed reinterpretation of this
candidate. Its g lapse residual on the unchanged rolling solution becomes
-4 identically in the normalized equations. Strict input controls reject
c=2, unadmitted branches, field
labels outside the interval, nonpositive roots, booleans and floats.

The independent Fraction engine directly differentiates the profile
quotients, evaluates full literal coframe polynomials and reconstructs
the lapse-margin coefficients from the primitive rational bounds. Its
finite fixtures check those algebraic identities, not the continuous
absence theorem in place of Sections 3–4.

## 6. What this says about the adopted V/B contract

The adopted V test starts with an actual vacuum and only then asks for
kinetic residues, non-tachyonic spectrum, amplitudes and dispersion
conditions. This action supplies no such proportional constant-clock
Minkowski vacuum within its stated local field domain. No positivity
amplitude or vacuum spectrum is claimed to have been computed here.

The contract explicitly permits a separately named smooth off-interval
extension. This proof neither assumes global real analyticity nor rules
out an extended action with a vacuum at another clock value. Extra
potentials, heavy-field backgrounds or changed coefficients are likewise
different candidates. A vacuum supplied by such a candidate would still
need the same-action, domain, source, scale and error controls of B.
The contract does not require a global physical-time trajectory or
scattering history from the bounce to the vacuum, and this proof imposes
neither. Original DHOST-row matching and original P8 remain open.
