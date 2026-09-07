# P8-S6.12.COMPOSITE — source-preserving symmetric vacuum-branch matching

This is a new bounded matching audit, not a revision of the original
DHOST classification or the frozen composite-bounce gates. It pins and
replays S6.11 (`34fa2ee35c2b33557e84b8fc181a202bed4793ed00d6aa8b76964ebbdc5dbe12`)
and the original CD/M1 witness and classification certificates.
The proportional GR/source results are prior art; see
[the normalization and attribution audit](notes/sources.md).

## Fixed gravitational model and two distinct matter applications

The quantitative vacuum model has constant positive Einstein coefficients
`G_parent=F_parent=M²`, `m4=M²m²`, `m>0`, `alpha=beta=1`, and bare
`beta_n=(0,0,1,0,0)`. The repository signature is `+---`. The physical
metric is always the prescribed composite metric, here called `mathcalG`:

`mathcalG=g(I+sqrt(g^-1 f))²`.

The exchange-parity identity extends to constant `beta_n=beta_(4-n)` with
the same equal Einstein coefficients and equal unit composite weights;
no arbitrary-beta vacuum health result is asserted.

Two matter applications are explicitly distinguished:

1. **Pinned one-scalar model.** The earlier composite family has one
   canonical scalar coupled to mathcalG, with a locally reconstructed
   potential. Its action and initial data are not CD/M1.
2. **Canonical-light M1 extension.** A separately named application of
   the symmetry theorem has canonical positive-field-metric clock fields
   psi^A coupled to mathcalG, and the separate exact free field
   `S_chi=integral sqrt|mathcalG| (partial chi)²/2`. Its potential and
   clock kinetic metric are independent of chi and there are no mixed
   chi kinetic couplings. Taking one clock gives the required two light
   scalars for a direct CD/M1 comparison. Adding this matter changes the
   earlier parent action/background problem; no pinned one-scalar
   solution is relabelled as a solution of this extension.

More generally the parity statement allows any prescribed light/source
action depending on mathcalG and light fields but not on the relative
coordinate. The reduced action retains that input exactly; it does not
erase higher operators if they were already inserted in the light action.
The zero-C/D conclusion applies to the specified minimally coupled,
first-derivative canonical-light class.

## Exact branch and conditional uniqueness

Use the regular positive-root Cayley chart

`Delta=(I-S)(I+S)^-1`, `S=sqrt(g^-1 f)`,

`g=mathcalG(I+Delta)²/4`, `f=mathcalG(I-Delta)²/4`.

Delta is mathcalG-self-adjoint. Required metric/root/chart invertibility,
signature and branch conditions are stated in the proof. They are not
inferred from a nonzero physical determinant alone.

The full action, with exchange-symmetric boundary terms or compactly
supported variations, is even in Delta. Thus Delta=0 is an exact
off-shell stationary relative-field branch, even with dynamical light
sources. On this branch the classical reduced action is

`integral sqrt|mathcalG|[-M² R[mathcalG]/4 -3M²m²/8 + L_light]`.

It induces no new C/D operators, curvature-squared terms or heavy-source
contacts at tree level. This is not a loop statement.

Calling this the unique eliminated branch requires specified zero-relative
initial/boundary data and a well-posed constrained problem, or a specified
stationary functional problem with the relevant linear operator an
isomorphism. The connected-branch uniqueness theorem is conditional on
these hypotheses. Positive algebraic mass or a potential Hessian is not
such an operator inverse. No uniform inverse on a rolling CD domain is
computed here. The conditional residual estimate `K*r/(1-eta)` does not
supply K or eta for that physical domain.

## Vacuum and direct matching verdicts

A proportional Minkowski candidate has `g=f=mathcalG/4`,
`M_eff²=M²/2`, and physical relative Fierz–Pauli mass `m_FP²=m²/4`.
For a fixed potential to realize this vacuum it must have
`V(phi_v)=-3M²m²/8`, `V_phi(phi_v)=0` and the specified nonnegative
canonical scalar mass Hessian. A strictly positive clock mass may be
chosen for the conditional vacuum test; a free chi remains massless.
These are point conditions on one action, not proof of an off-tube
extension matching the pinned rolling reconstruction. The constant
vacuum term shifts all five beta_n by `-3/8`; using bare beta2 in a
standard mass formula gives the wrong mass.

For the same prescribed physical metric and clock/operator dictionary,
the exact zero-relative canonical-light branch cannot reproduce the old
CD/M1 covariant action. In pinned `M=tau=1` normalization the target has
`F2_X(0,1)=-1/2`, `A3(0,1)=1`, while the branch has both zero.
A remainder claimed to repair those jets must therefore have respective
absolute magnitudes at least `1/2` and `1` there.

Independently of that operator dictionary, a spatially flat homogeneous
solution of the branch with canonical NEC matter has nonincreasing
physical H. It cannot approximate the old CD H at both
`T=+-tau/2` with endpoint errors strictly below `8/(5tau)`.
For the exact CD background and exact free M1 charge, the null-equation
residual needed to repair the GR branch has magnitude at least
`(801/100)M_eff²/tau²` at the bounce. Approximate acceleration/chi-speed
budgets are given in the proof. These are necessary discrepancies,
not computed quantum or heavy-state corrections.

No full-parent, all-state, all-bimetric, loop, cutoff or UV exclusion is
claimed. Nonzero relative condensates, different boundary states, loss
of functional invertibility, symmetry breaking and other parent actions
remain separate possibilities. S6 does not require a physical-time
trajectory from vacuum to bounce; none is assumed. The original row
classification remains intact and S6/P8 remain open.
