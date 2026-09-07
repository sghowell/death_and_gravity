# P8 S6.18.DETERMINANT — regular determinant-sum backgrounds cannot bounce

This is an exact background exclusion for a named action and domain, not a
UV verdict, stability result, general Lorentz-branch theorem or closure of
the original C/D matching problem. It does not edit the frozen STAR result.

## 1. Action, sources and regular domain

Let a finite nonempty set of metrics have positive constant Einstein
coefficients `G_i`, constant real cosmological terms `Lambda_i`, and the
rank-one interaction

```
S = -sum_i G_i/2 int sqrt|g_i| (R_B[g_i]+2 Lambda_i)
    -lambda int det U + sum_i S_m,i[g_i,psi_i],
U=sum_i beta_i e_i,       lambda=m^4>0.
```

The real beta coefficients are constant and need not have the same sign.
Each matter sector couples minimally only to its own actual EH metric;
the sectors do not share fields or directly couple to one another. On the
background their stresses are homogeneous/isotropic, separately conserved,
and have `nu_i=rho_m,i+p_m,i>=0`. Positive-metric canonical scalar fields
are a sufficient example. No NEC is imposed on the interaction stresses.

Use the P8(b) convention `+---`, `R_B=-6(DH+2H^2)`, so
`3G_i H_i^2=G_i Lambda_i+rho_m,i+rho_int,i` in spatially flat FLRW. The
null Einstein equation is `-2G_i D_i H_i=nu_i+rho_int,i+p_int,i`.

The coframes are simultaneously aligned and diagonal on a connected open
flat-FLRW chart: `e_i=diag(n_i,a_i,a_i,a_i)`, with all actual `n_i,a_i`
smooth, finite and strictly positive. A common Lorentz transformation is
irrelevant; no uniqueness claim about other relative-Lorentz branches is
made. Define

```
S_a=sum_i beta_i a_i,        S_N=sum_i beta_i n_i.
```

The regular-sum theorem assumes `S_a*S_N!=0` throughout the interval. All
nonzero beta having the same sign is sufficient for this with positive
coframes; it is not necessary for the algebraic proof. That proof does not
establish the physical health of mixed-sign parameter choices.

## 2. Physical result

Choose any actual matter metric `g_r`. If `beta_r!=0`, sum only over the
active component `I={i:beta_i!=0}` and use its physical proper time:

```
dT=n_r dt, H=D_T log a_r, y_i=a_i/a_r, c_i=n_i/(n_r y_i),
y_r=c_r=1, K=sum_I G_i y_i^2>0.
```

The full background equations imply, including at zeros of H,

```
K H'-H K'/2=-sum_I c_i y_i^4 nu_i/2,
(H/sqrt K)'=-sum_I c_i y_i^4 nu_i/(2 K^(3/2))<=0.
```

Consequently `H(T0)<=0` implies `H(T)<=0` at every later time in the same
regular interval. Physical contraction-to-expansion, including a degenerate
transition, is excluded. No principal-cone bound, vacuum spectrum, auxiliary
TT inverse or expansion in a source parameter enters this argument.

If `beta_r=0`, that metric is disconnected and independently obeys the
usual flat GR null equation `H'=-nu_r/(2G_r)<=0`. It is not inserted into
the coupled K and does not need regularity of an unrelated sum. All-zero
beta is simply a collection of separately coupled GR metrics.

## 3. Matter frame and remaining boundary

On an invertible stationary auxiliary branch,
`L_aux=lambda det w [3-tr(w^-1 U)]` gives `w=U` and `-lambda det U`.
This is the old auxiliary potential with `B=-3lambda/2` and
`p_i=lambda beta_i/2`. The matter actions stay on the original EH leaves,
not on w. Substituting a new matter action on w changes its equation and
does not retain this elimination. A positive w coframe chart exists only
where the signs of the summed entries permit it; it is not assumed for the
direct determinant proof.

The polynomial interaction is evaluable at singular U, but the regular
proof must not divide by `S_N=0` there. A recorded Bianchi-only point with
`S_N=0,S_a!=0` fails velocity locking. It is explicitly **not** an actual
solution, viable branch or bounce. No verdict is given on persistent
singular-sum strata, their rank/constraint health, or a continuation across
them. Spatial curvature, nonaligned Lorentz branches, variable coefficients,
derivative mixing, direct composite/multimetric matter and physical NEC
violation also remain outside this theorem. No high-frequency-cone result
is substituted for a physical light-EFT or full UV verdict.

The exact symbolic and independent Fraction/covariant tests support the
written proof, not an all-orders proof-assistant formalization. S6.16,
S6.13 and the adopted S6 contract are pinned and replayed read-only; no
concurrent TREE certificate is used.
