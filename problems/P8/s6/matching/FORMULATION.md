# S6.2 — First parent-EFT matching tests

Opened 2026-09-04 under the adopted S6 contract and authorization to proceed.
The following are selected **matching ansatzes**, not additional universal
UV assumptions. Their analysis followed exploratory algebra; this is not a
preregistered blind search. No old witness or certificate is changed.

## A. Stable algebraic heavy-scalar matching

Start with the positive-lambda vacuum of S6.1 and its massive-mediator
regression. Test finite collections of real heavy scalars h^a, whose local
zero-heavy-derivative Lagrangian, at fixed physical clock phi, is

```
L_alg = Z(phi,h)*X/2 - V(phi,h),  X=(partial phi)^2.
U=V-Z*X/2; U_a(phi,h_*(phi,X),X)=0
K_ab=U_ab|_* positive definite.
```

This is an affine-in-X parent before elimination, with a stable stationary
minimum. The positivity of K is a stipulated property of this branch. It
must not be inferred just from a nonzero physical heavy-mode frequency in
a mixing background. General kinetic mixing and retained heavy derivatives
must be included in a remainder rather than discarded without an estimate.

In the fixed P8 physical metric and clock variable, test
`F_target=P_alg+R`, where `P_alg=L_alg(phi,h_*(phi,X),X)` and R includes all
omitted effects that feed this coefficient. Compare the **second X jet**,
not only function values. No derivative field redefinition or integration-
by-parts dictionary is assumed to preserve F_XX; if one is used, it needs
an explicit matching dictionary and error estimates.

A curvature source linear in R may be included in U and frozen to R=0 for
this particular coefficient test, but that extension alone does not prove
a metric/curvature matching result. No loop or general UV exclusion follows
from failure of this tree-level algebraic ansatz.

## B. Minimal classical parent, checked without derivative expansion

Test the full classical action

```
S_parent=integral sqrt(-g)*[-M^2 R/2
             +G_AB(Phi)*partial Phi^A.partial Phi^B/2 - V(Phi)]
G_AB positive definite; M>0 constant.
```

The metric g is the physical P8 metric and spatial sections are flat. All
background fields are homogeneous. No quantum expectation stress, additional
NEC-violating source, or higher-curvature operator is inserted. The kinetic
metric can be curved and mixed, and the number of fields/potential is arbitrary.
This is a test of an actual parent solution, not a claim that every healthy
UV theory belongs to this class.

## C. Nonminimal scalar-curvature parent

Also test `-Q(Phi) R/2`, Q>0, in place of the constant Einstein coefficient.
Q is a function of the parent fields, not of their derivatives. Require the
Einstein-frame field metric to be positive definite; the Jordan-frame field
metric itself need not be positive. Keep the original metric g physical.
The auxiliary transformation `g_E=(Q/M^2)*g` is used only to diagnose the
parent equations, not to identify different physical-matter-frame witnesses.

The target D tube has Q_target=M^2. For a parent intended to match this at
the bounce t=0, put q=Q/M^2 and bound its **physical-time** jets by

```
|q(0)-1|<=epsilon0<1,
|tau*q_dot(0)|<=epsilon1,
|tau^2*q_ddot(0)|<=epsilon2.
```

Derive a necessary relation among these errors rather than postulating that
small amplitude error implies derivative control. These bounds concern Q
evaluated on the parent background; a claimed EFT matching must relate it
to the effective Planck coefficient explicitly. An order-one discrepancy,
rapid heavy motion or a different field/operator map is not a small error
in this ansatz. The exact target scale factor is a=1+(t/tau)^2; an optional
extension bounds an error in H_dot itself separately.

## Evidence and outcomes

Use exact parent variation, an independent conformal-geometry derivation,
stationary-elimination identities, the pinned D witness and negative controls.
Generic positive-matrix/calculus implications are written proofs, not
formalized by a few symbolic matrix examples. Known parent NEC and Schur
arguments are not claimed as new discoveries.

Record whether these ansatzes can match the chosen tube and how large any
necessary correction is. Failure is not a D-row exclusion and cannot be
turned into a positivity verdict on the rolling DHOST scalar. More general
derivative/curvature parents, other witnesses, nonadiabatic branches and
quantum effects remain separate questions.

## D. Regular metric maps and the two tails

As a follow-up to C, test whether a large but regular change of metric can
repair the matching. For an Einstein/positive-field-metric parent, retain its
exact two-derivative tensor action and Einstein-frame NEC. A homogeneous
disformal map has `a_E=sqrt(C)*a`, `dT_E=sqrt(W)*dt`, C,W>0 at every finite
time. The fields are homogeneous, so no extra tensor principal term comes
from their minimal two-derivative kinetic action. Match physical tensor
coefficients, not just a bare curvature coefficient. Permit arbitrarily
large finite derivatives of C and W. Test both the required nonzero tensor
tails and a finite-slice tensor-value error budget at t=(-L*tau,0,L*tau).
The target scale factor remains exact in this finite-slice test.

Separately test exact, open-tube changes of variables of the D effective
action itself to a **quartic Horndeski** action, not necessarily an Einstein
or NEC parent. Use smooth one-clock maps
`g_E=C(phi,X)*g+B(phi,X)*dphi*dphi`, `W=C+B*X`, with C,W positive and
`J_map=C-X*C_X-X^2*B_X` nonzero wherever a map is claimed. These functions
are unrelated to the optional function-group names C and D. Derive necessary
conditions from lapse-velocity mixing and the seed Horndeski tensor identity.
Do not identify a tensor-normalized frame with an Einstein/conventional-
scalar action, or confuse a nondegenerate metric with an invertible map.

This second test concerns an exact field map, not integrating out heavy
fields, an approximate change of action, all higher-order Horndeski theories,
or arbitrary derivative redefinitions. A failed auxiliary map does not
invalidate the already-certified original D variables. Both tests are
selected matching routes, not new universal UV assumptions.
