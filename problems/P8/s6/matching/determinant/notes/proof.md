# Literal determinant-sum no-bounce proof

The theorem uses exactly the [formulation](../FORMULATION.md). A prime is
the proper-time derivative of the chosen actual metric `g_r`; dots initially
mean the common coordinate derivative. The null density of actual matter
is denoted `nu_i` to distinguish it from the lapse `n_i`.

## 1. Action variation before fixing a lapse

On aligned flat-FLRW coframes, the full non-pairwise interaction is
`L_int=-lambda S_N S_a^3`. Vary every lapse and scale independently before
fixing `n_r=1`. In the actual metric i,

```
rho_int,i=-a_i^-3 partial_(n_i) L_int
         =lambda beta_i S_a^3/a_i^3,
p_int,i=(3 n_i a_i^2)^-1 partial_(a_i) L_int
       =-lambda beta_i S_N S_a^2/(n_i a_i^2).             (1)
```

The independent engine also varies all sixteen entries of each coframe
before specializing to diagonal isotropy. Each diagonal spatial variation
gives (1), and off-diagonal first variations vanish at the aligned point.
The full Lorentz equations are consequently satisfied on this branch;
other branches are not classified.

The individual interaction null densities can have either sign, but they
obey a literal polynomial cancellation without division by either sum:

```
sum_i n_i a_i^3 (rho_int,i+p_int,i)
 =lambda S_a^2 [S_a sum_i beta_i n_i-S_N sum_i beta_i a_i]
 =0.                                                     (2)
```

Endpoint cosmological terms have zero null density. Einstein variation in
the declared B convention gives
`-2G_i D_i H_i=nu_i+rho_int,i+p_int,i`, where
`D_i=n_i^-1 d/dt`, `H_i=dot a_i/(n_i a_i)`.

Because actual matter is separately conserved and every `G_i` is constant,
each Einstein equation implies conservation of its interaction stress.
Substituting (1), without a Hubble division, gives

```
D_i rho_int,i+3H_i(rho_int,i+p_int,i)
 =3lambda beta_i S_a^2/(n_i a_i^3)
    [dot S_a-S_N dot a_i/n_i]=0.                          (3)
```

If `beta_i=0`, both interaction stresses and (3) vanish identically; that
metric must be omitted from locking. For nonzero beta, `lambda>0`, positive
coframes and `S_a S_N!=0`, (3) implies

```
dot a_i/n_i=dot S_a/S_N=Q(t)   for every active i.         (4)
```

Equation (4) is valid at a simultaneous zero of every scale-factor
velocity. The coordinate-lapse quotient `dot a_i/dot a_r` is unnecessary.
The proof never divides by H, Q or an actual matter density.

## 2. Physical proper clock, positive kinetic weight and cancellation

Fix active r and put `dT=n_r dt`, `a=a_r`, `H=a'/a`,
`y_i=a_i/a`, `c_i=n_i/(n_r y_i)>0`. In particular `y_r=c_r=1`.
The common velocity in (4) gives the exact kinematic formulas

```
H_i=H/y_i,
y_i'=y_i(c_i-1)H,
D_i=(c_i y_i)^-1 D_T,
D_i H_i=[H'-(c_i-1)H^2]/(c_i y_i^2).                    (5)
```

Multiply the null equation for i by `c_i y_i^4`, not by `y_i^2` or an
unweighted unit. In terms of the coordinate variables this factor is
`n_i a_i^3/(n_r a^3)`. Hence (2) cancels every interaction term, leaving

```
sum_I G_i y_i^2 [H'+(1-c_i)H^2]
          =-sum_I c_i y_i^4 nu_i/2.                       (6)
```

There are finitely many active indices; their Einstein coefficients and
ratios are positive. Set

```
K=sum_I G_i y_i^2>0,
K'=sum_I 2G_i y_i y_i'=2H sum_I G_i y_i^2(c_i-1).
```

This is an actual chain rule, with K' allowed either sign. Replacing the
sum in (6) yields

```
K H'-H K'/2=-sum_I c_i y_i^4 nu_i/2,
(H/sqrt K)'=-sum_I c_i y_i^4 nu_i/(2K^(3/2))<=0.          (7)
```

This exact identity applies on the entire stated regular interval. Integrate
it between any two interior times. The sign of H equals that of
`H/sqrt K`, so a nonpositive H cannot later become positive. This covers
degenerate transitions too. At H=0 the undivided equality in (7) gives
`K H'=-sum_I c_i y_i^4 nu_i/2<=0`, but the integrated identity, not that
pointwise observation alone, proves no crossing.

If all actual sources except r vanish, the numerator reduces to `nu_r`.
If several metrics have matter, their positive clock/measure weights remain;
no artificial single-source identification is made. Only their separate
conservation and NECs are used. No null condition on an effective spin-2
fluid or subluminality of a principal metric cone is assumed.

For `beta_r=0`, the r interaction vanishes identically and its own Einstein
null equation gives flat GR monotonicity directly. Including disconnected
fields in the active K would incorrectly import (4) into equations `0=0`.

## 3. Actual controls, including a rolling source

To show the regular action/domain is nonempty, take three actual metrics,
`lambda=1`, `beta=(1,-1,1)`, `G=(1,2,3)`, and constant ratios
`y=(1,2,3)` in the physical clock of the first metric. Then `S_y=2`,
and the interaction densities are `(8,-1,8/27)`. Define

```
Lambda=(-8,1/2,-8/81),
a_i=y_i T^(1/3), n_i=y_i,
phi_i=sqrt(2G_i/3) log(T/T_*), V_i=0, T>0.
```

Here fixed units set the harmless reference scales to one; `T_*>0` simply
makes the logarithm dimensionless. The actual scalar densities and
pressures are `rho_m,i=p_m,i=G_i/(3y_i^2 T^2)`, and the exact scalar
current equation `(a_i^3 phi_i'/n_i)'=0` holds. Each interaction has
`p_int,i=-rho_int,i`; the chosen constant Lambda cancels its density in
that metric's Einstein equation. Both independent Einstein equations for
all three metrics therefore hold. The actual summed coframe is positive
and invertible: `S_N=2`, `S_a=2T^(1/3)`.

The physical kinetic weight is K=36, giving
`(H/sqrt K)'=-1/(18T^2)`. This is an expanding exact example, not a bounce,
not original C/D matching, and not a perturbative-health claim about mixed
signs. It is a test of the broader algebraic theorem, not of the paper's
same-sign restriction.

A second family has the same beta,G,y and
`a_i=y_i exp(H0 T)`, `n_i=y_i`, with zero real matter and
`Lambda_i=3H0^2/y_i^2-rho_int,i/G_i`.
All Einstein equations again hold. H0 is any fixed real number, and the
normalized derivative vanishes. This covers expanding/contracting flat
de Sitter patches and the H0=0 flat case, with no assertion about their
global extensions. It checks the allowed NEC-saturation limit.

## 4. An exact auxiliary identity, with the source not moved

For invertible w consider the polynomially extendible expression

```
L_aux=lambda det w [3-tr(w^-1 U)]
     =-2 det w [B+tr(w^-1 sum_i p_i e_i)],
B=-3lambda/2,  p_i=lambda beta_i/2.                       (8)
```

Vary w freely as a matrix. After multiplication by the invertible matrix
factors its Euler equation is
`(3-tr Q)I+Q=0`, where `Q=w^-1 U` (or its transpose according to the
matrix-gradient convention). Taking the trace gives `tr Q=4`, hence
`Q=I`, `w=U`; the stationary action is `-lambda det U`.
The independent coframe/matrix audit checks the full stationary variation,
not only a trace fit on one background.

The original matter actions depend only on their EH leaves and not on w,
so (8) is source-preserving for the present theory. By contrast, adding
`S_m[w,psi]` sources the w equation and generally changes its solution.
It would be a different composite/auxiliary-matter model; the current paper
and this theorem do not automatically apply to it. In a positive w chart
one additionally needs the summed coframe entries to admit that chart.
The direct proof (1)–(7) uses only positivity of actual metric lapses/scales
and nonzero sums, not a positive auxiliary chart.

## 5. Sum degeneracies and precise non-conclusions

Equation (3) was retained before division. If `S_N=0` persists while
`S_a!=0`, it yields only `dot S_a=0` for active indices. Pressures in (1)
are then zero but densities can be nonzero and signed. Velocity locking
does not follow. For example the point

```
beta=(1,1,-1), a=(1,2,1), n=(1,1,2), dot a=(1,-1,0)
```

has `S_a=2`, `S_N=dot S_a=0`, and every unfactored Bianchi residual zero,
but `dot a_i/n_i=(1,-1,0)` are distinct. This is a **Bianchi-only control**:
no Einstein/source equations or existence/health of a persistent branch
are asserted. It prevents a false extension of (4) by algebraic division.

At `S_a=0`, (1) vanishes; this alone must not be used to dispose of the
different `S_N=0,S_a!=0` stratum. The determinant polynomial exists at
singular sums, but the present regular theorem does not classify these
branches or give a matching/UV verdict on them. The source paper assumes
an invertible sum and common aligned Lorentz branch; the statement here
retains those relevant restrictions instead of claiming they follow from
general physical consistency.

The finite exact tests verify the action, clock, source weights and stated
controls. The mathematical sign conclusion is the written integration of
(7), with explicit open-domain hypotheses. It is not a statement about
arbitrary interacting spin-2 graphs, nonaligned Lorentz branches, curved
spatial slices, variable couplings, quantum NEC violation or changed matter
couplings. Nor does a background result establish a tensor/scalar/vector
cutoff, a high-frequency light-mode matching, or any UV completion.
