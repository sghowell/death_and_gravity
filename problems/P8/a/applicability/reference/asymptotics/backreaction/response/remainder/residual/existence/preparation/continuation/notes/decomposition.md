# The complete A.11 map, with one block resummed exactly

The distinction is essential: the block below is not the full Frechet
linearization at the actual time-dependent metric/state. It can be
isolated **without dropping a term** from the original fixed-point
equation; the remaining functional must still be estimated.

## 1. Exact trace identity

Use u=-a''/a, h=a'/a, q=S/a^2, h'=-u-h^2 and

    P[h,u]=u^2/4+h^2(u+h^2)/30.

Then

    a^2(q''+2hq')=S''-2hS'+2(u+h^2)S,
    (a^2 q')'=a^2u/(60delta)-P[h,u]

on the source-free interval. The two expressions reproduce the entire
A.11 trace, including the anomaly and the positive Einstein coefficient.
Ordinary radiation is traceless but its fixed normalization remains in
the initial density constraint; a trace calculation cannot replace that
constraint.

## 2. Shared-history variation and auxiliary variables

Translate the start of a proposed continuation interval to t=0. Choose
a smooth baseline continuation ub,ab of the already known past. This
may be an approximation, not an exact solution. Do not reset modes at
t=0. For any candidate, set W=u-ub=I X, X=u'-ub'; the past difference is
zero. Common metric initial data determine a through a''=-ua. The actual
quantum modes and R[u] are evolved from the original common free past.

Solve the **auxiliary** integrated trace for q on each candidate metric,
with q(0),q'(0) equal to the actual initial Wick data S/a^2 and its
derivative. Define qb in precisely the same way on the baseline. Unless
it is a fixed point, this auxiliary q is not asserted to equal S[u,a]/a^2.
This separation is how A.11 avoids reintroducing an uncontrolled derivative
of the singular logarithmic operator into its Banach map.

Let d_f=d(af) and write

    S[u,a]=Lf u+N[u,a],
    N[u,a]=R[u]+(d(a)-d_f)u,
    Lf=4pi^2 D_gammaE+d_f Id.

On variations supported after t=0 the singular operator has no new past
endpoint term. Since the past histories coincide, subtracting their
original retarded integrals is exactly the translated causal operator
on W. The nonlinear difference still knows the original state history.
For smooth flat variations, all derivative commutations are literal;
the same resummed linear block is defined causally for continuous input,
without claiming an automatic C1 output.

## 3. Every term of the resummed map

Set

    r_b=(ab^2 qb)'-S[ub,ab]',
    Delta N'=R[u]'-R[ub]'+(d(a)-d_f)X
             +(d(a)-d(ab))ub'-(hu-hb ub)/2.

Here r_b is the actual baseline fixed-point residual, not set to zero.
The auxiliary identity yields

    Delta[(a^2q)']
      =I[(a^2u-ab^2ub)/(60delta)-Delta P]+2 Delta(a^2 hq).

Consequently equality of the actual and auxiliary Wick derivatives is
**exactly**

    (Lf-c I^2/2)X=Gf[X], c=af^2/(30delta),
    Gf[X]=r_b
       +I[(a^2-af^2)u-(ab^2-af^2)ub]/(60delta)
       -I[P[h,u]-P[hb,ub]]+2(a^2hq-ab^2hbqb)-Delta N'.

Common actual initial Wick data then give equality of the Wick functions,
not merely their derivatives. The density constraint must still hold and
propagate by the full conserved stress, as in A.11. If a preparation
source is kept on, its exact conserved-source terms must also remain in
Gf; the present decomposition is deliberately for the unforced future.

This formula isolates both the local finite coefficient and the whole
constant Einstein term. The nonlinear actual-state derivative, rolling
Einstein mismatch, anomaly and q geometry are not discarded. Replacing
Gf by a prescribed function defines the linear comparator analyzed next.
Setting the remaining terms to zero is **not** a theorem that the rolling
state has the flat-space Fréchet derivative. In particular af constant,
h=u=0 with the original nonzero ordinary radiation is not a full static
solution of the density constraint, even though its linearized trace
contains this block.
