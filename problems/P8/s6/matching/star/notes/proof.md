# Proof of the constant-star background exclusion

All primes below mean the physical proper-time derivative `D_T`, except
derivatives of the displayed polynomials, which are explicitly `d/dR`.
The hypotheses are exactly [FORMULATION.md](../FORMULATION.md).

## 1. Literal interaction and the two proper clocks

For one leaf suppress its index. The square-root eigenvalues are
`(N,R,R,R)` with `N=n_u/n_i`, `R=a_u/a_i`. Expanding their elementary
polynomials before variation gives

```
A(R)=beta0+3 beta1 R+3 beta2 R^2+beta3 R^3,
B(R)=beta1+3 beta2 R+3 beta3 R^2+beta4 R^3,
J(R)=beta1+2 beta2 R+beta3 R^2,
sum beta_n e_n(N,R,R,R)=A(R)+N B(R),
L_int=-2 a_i^3 [n_i A(a_u/a_i)+n_u B(a_u/a_i)].
```

The lapse and scale variations use, for either metric,
`rho=-a^-3 partial_n L_int` and
`p=(3 n a^2)^-1 partial_a L_int`. Hold the other lapse/scale fixed; do not
fix both lapses before varying. Using `dA/dR=3J` and `3B-R dB/dR=3J`,

```
rho_i=2A,              p_i=-2[A+(N-R)J],
rho_u,i=2B/R^3,        p_u,i=-2J/(N R^2)-2(dB/dR)/(3R^2),
n_i,int=2(R-N)J,       n_u,int=2J(1-c)/R^3,       c=R/N,
n_u,int+n_i,int/(N R^3)=0.                              (1)
```

Here `n_i,int` denotes an interaction null density and is not the lapse
`n_i`. In particular no NEC is imposed on either interaction stress.
The relative signs in (1) are compatible with the known reciprocal
interaction-NEC relation, but here they follow from the actual action.
Independent tests vary all three spatial coframes separately before
specializing to isotropy; no anisotropic equation is silently discarded.

Removing the usual Einstein boundary term gives the flat minisuperspace
Lagrangian `-3G a dot(a)^2/n`. Its lapse and scale equations yield
`3G H^2=rho` and `-2G D H=rho+p`. Thus, for matter only on the center,

```
3G_i H_i^2=2A_i,
G_i D_tau_i H_i=(N_i-R_i)J_i,
G_u H_u'=-n_h/2-sum_i J_i(1-c_i)/R_i^3.                 (2)
```

These formulas also hold for `G_u=0` as an algebraic central null equation.
A separate central cosmological term cancels in the last equation.

Kinematically, without any equations of motion,

```
R_i'=R_i H_u-c_i H_i,
D_tau_i=(R_i/c_i)D_T,
D_tau_i R_i=R_i(N_i H_u-H_i).                           (3)
```

The leaf has no independent matter, so its full Einstein equation implies
conservation of its interaction stress. Substituting (1) and (3) gives the
unfactored identity

```
D_tau_i rho_i+3H_i(rho_i+p_i)
       =6 N_i J_i(R_i)(R_i H_u-H_i)=0.                  (4)
```

No division by `J_i`, `H_u`, a null density or a TT Hessian is involved.
Only strictly positive regular lapses and ratios will be divided below.

## 2. A fixed dynamic stratum

When `J_i(R_i)!=0`, (4) implies `H_i=R_i H_u`. Differentiate this identity
on its open domain and apply (2),(3):

```
R_i'=R_i(1-c_i)H_u,
(G_i/R_i^2)[H_u'+(1-c_i)H_u^2]=J_i(1-c_i)/R_i^3.       (5)
```

On a region with a fixed set D of dynamic leaves and all other genuine
leaves identically at roots, the latter interaction null densities vanish.
Adding (5) for D to the central null equation therefore gives

```
K_D=G_u+sum_D G_i/R_i^2,
K_D'=-2H_u sum_D (G_i/R_i^2)(1-c_i),
K_D H_u'-H_u K_D'/2=-n_h/2.                            (6)
```

Whenever `K_D>0`, division by `K_D^(3/2)` gives the stated nonincreasing
`H_u/sqrt(K_D)`. At `H_u=0`, (6) remains an undivided equation
`K_D H_u'=-n_h/2<=0`. No restriction on a link sign or on `c_i-1` was used.
Equation (6) is not asserted with algebraic leaves added to D, or with
`K_D=0`.

## 3. Algebraic interiors and arbitrary zero-set boundaries

For a genuine link the constant polynomial `J_i` is nonzero and has at most
two roots. On a connected open interval where `J_i(R_i(T))=0`, continuity
and the finite root set force one constant positive `R_i`. Its leaf null
equation in (2) gives `D_tau_i H_i=0`, hence constant `H_i`, and (3) gives

```
H_u=c_i H_i/R_i,       H_u'=(c_i'/c_i)H_u.              (7)
```

The sign cannot change on such an interval; a zero constant `H_i` gives
`H_u=0`. This is an actual background consequence, not a claim about scalar
or vector perturbations on that branch.

We next avoid any tacit assumption that there are only finitely many
switches. Let `Z_i={T:J_i(R_i(T))=0}`, a closed subset of the regular open
time interval. Let `F_i=H_i-R_i H_u`, which is continuously differentiable.
On the open complement of `Z_i`, (4) gives `F_i=0` and then `F_i'=0`.
Every boundary point of `Z_i` is a limit of this complement. Continuity of
both `F_i` and `F_i'` gives their vanishing at that boundary point too.
Consequently (5), including its differentiated relation, is valid at every
point outside `int Z_i`. Equation (7) holds locally at every point of
`int Z_i`. This covers isolated roots, accumulating roots, and zero sets
with empty interior without an auxiliary rank assumption.

If `J_i` is identically the zero polynomial, the link has only beta0 and
beta4 endpoints. These are separate cosmological terms, not an interaction
link. Omit that leaf in the following argument; no spurious Bianchi
constraint is inferred from `0=0`.

## 4. A compact positive-part argument through all strata

Fix an arbitrary compact physical-time interval `[T0,T1]` lying in the
regular open chart. Positivity and smoothness make the following constant
finite, with the maximum over the finite set of genuine links:

```
C=max_i {sup_[T0,T1] |R_i'/R_i|, sup_[T0,T1] |c_i'/c_i|}.
```

Take `C=0` for an empty genuine set. We will prove
`H_u'<=C H_u` at every point where `H_u>=0`.

If some genuine link has that point in `int Z_i`, its local (7) implies
this inequality directly, irrespective of all other branch assignments.
If no genuine link has an algebraic interior there, section 3 establishes
(5) at the point for every genuine link. Use
`K=G_u+sum_genuine G_i/R_i^2`, which is strictly positive by the formulation.
The pointwise combination in (6) is therefore valid there. Independently
differentiate this smooth positive K by the chain rule:

```
K'/(2K)=-sum_i (G_i/R_i^2)(R_i'/R_i)/K.
```

The right side has absolute value at most C because the weights are
positive and their sum is at most K; the remaining central weight has zero
rate. Thus (6), `n_h>=0`, and `H_u>=0` again give `H_u'<=C H_u`.
If there are no genuine links and `G_u>0`, this same argument is just
`G_u H_u'=-n_h/2`.

Set `Y=max(H_u,0)`. It is absolutely continuous. On `H_u>0`, its almost
everywhere derivative is `H_u'<=C Y`; on `H_u<0` it is zero. The standard
chain rule for the positive part gives zero derivative almost everywhere
on the zero level set as well. Therefore `Y'<=C Y` almost everywhere.
The absolutely continuous integrating factor obeys
`(exp[-C(T-T0)]Y)'<=0`. If `H_u(T0)<=0`, nonnegativity and `Y(T0)=0`
force `Y=0` throughout `[T0,T1]`. Exhausting compact future subintervals
proves the theorem. This proves the sign result even for a degenerate
transition; merely checking `H_u'<=0` at isolated zeros would not suffice.

The coefficient C need not be uniformly bounded at the edge of the open
chart. It exists on each compact regular interval. We neither cross a
singular endpoint nor infer completeness/incompleteness from this theorem.
The argument never uses a globally exact K through an algebraic interior.

## 5. Two actual controls and the original parent map

An explicit nonempty algebraic solution is, for `T>0`,

```
G_u=G_i=1,  beta=(0,1,-1/2,0,1/2),
a_u=a_i=T^(1/3),  n_u=1,  n_i=1/(3T),
phi=sqrt(2/3) log T,  V=0.
```

It has `R=1`, `A=3/2`, `B=J=0`, `H_i=1`, `H_u=1/(3T)`.
The actual leaf stress is `(rho_i,p_i)=(3,-3)`; the central interaction
stress vanishes. The canonical scalar has `rho_m=p_m=1/(3T^2)` and
`(a_u^3 phi')'=0`. Both independent Einstein equations for both metrics
hold. The correct dynamic coefficient is only `K_D=G_u=1`.
Using the wrong all-leg `K=2` gives
`2K H_u'-H_u K'+n_h=-2/(3T^2)`, nonzero. Equivalently the error in the
normalized monotonicity formula is `-1/(6 sqrt(2) T^2)`. This control is
expanding and does not contradict the sign theorem.

For the excluded disconnected center take all beta zero, `G_u=0`, one
`G_i=1` Minkowski leaf and zero scalar. The interaction action and every
one of its variations vanish. The leaf Einstein and scalar equations
vanish, while the center has no kinetic equation. Thus
`n_u=1`, `a_u=1+T^2` is genuinely undetermined and `H_u'(0)=2>0`.
The replay derives the zero potential jets and Einstein/source residuals;
it does not merely accept displayed zero constants. Adding `G_u>0`
would impose `3G_u H_u^2=0` and destroys this control.

For the frozen original auxiliary model use two leaves with
`beta_i3=p_i`, `beta_i0=b_i`, all other leaf betas zero, `G_u=0`, and
the separate central potential `-2B_central sqrt|h|`. Then
`J_i=p_i R_i^2` and `rho_u,i=6p_i/R_i` reproduce the literal S6.13 action.
Every nonzero p is genuine and has no positive algebraic root. At least
one such link suffices for the sign obstruction. The theorem is broader
than that original model but is still a named action-class exclusion.

## 6. Verification and scientific scope

The symbolic engine varies the actual homogeneous action and checks its
identities. The independent engine uses Fraction arithmetic, polynomial
coefficients, first-order jets of four independent coframe entries and
proper-time fixtures; it imports no primary p8_star formulas. The reserved
covariant audit independently verifies the same action, fixture equations
and compact comparison identities. Domain tests reject binary floats,
nonfinite values, wrong signs and missing kinetic coefficients. Corrupted
certificate controls test source and claim integrity.

These finite exact replays support the written all-interval proof; neither
pytest nor polynomial identities alone encode its topology and Gronwall
steps. No assertion about perturbative stability, a heavy cutoff, a
scattering prescription or finite-error matching is inserted. A different
action or non-NEC source must be analyzed separately under the adopted S6
contract. There is no universal DHOST-row or UV-completion verdict here.
