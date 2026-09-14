# Conditional fixed-source Ward and Hessian transport

This is a conditional identity for a separately specified,
differentiable regulated functional. It does not assert that the
complete P8 continuum functional already satisfies its hypotheses.

Use Minkowski conventions
Z_Psi[J]=integral mu rho exp(i(S+sPsi+J_a F^a)/hbar),
W=-i hbar log Z, phi^a=delta W/delta J_a, and
Gamma=W-J_a phi^a. The original observables F^a are held fixed
before reduction; they can be nonlinear and composite. No new
linear source in a different field coordinate is substituted.

Let A=deltaPsi be odd. The graded integration identity is

< sA > = (i/hbar) J_a < A sF^a > + E[A].

Here E[A] contains the variation of the measure, regulator, state,
counterterms or boundary density, plus integration-boundary fluxes.
More explicitly, if the Lie variation of the source-free density
w is A_defect w, then integrating the graded derivative of
A w exp(i JF/hbar) gives
E[A]=< A A_defect > plus the normalized boundary term.
The sign follows from the minus in the odd Leibniz rule.
Explicit variations of non-gauge-fermion action or state data,
if present in a proposed construction, must also be included.

Thus, for fixed original sources and a variation only of Psi,

deltaW=J_a K^a+E,
K^a=(i/hbar)<deltaPsi sF^a>,
deltaGamma+Gamma_,a K^a=E.

At vanishing defects a stationary effective-action VALUE is
gauge independent. The coordinate mean is not generally fixed.
Differentiating the stationary equation yields delta phi=K only
where the relevant stationary Hessian complement is invertible;
otherwise only the corresponding Hessian equation is justified.
With nonzero E its source/field derivatives remain.

Twice differentiating the identity gives

delta Gamma_,ij + K^a Gamma_,aij
 + K^a_,i Gamma_,aj + K^a_,j Gamma_,ai
 + K^a_,ij Gamma_,a = E_,ij.

The last onepoint contact on the left is generally nonzero off
shell. On a stationary mean path with E=0 it drops out, leaving
dH=-K_,phi^T H-H K_,phi. A regular finite transformation therefore
acts by congruence and preserves inertia. This does not establish
regularity of the actual quantum gauge flow, nor allow dropping
a onepoint of only part of the action. The P8 nonreference
principal obstruction and unknown physical cutoff are unchanged.

## Exactly integrated control

For an independent finite example let q be a real Gaussian of
variance v>0, let y be a gauge coordinate, and choose
chi=y-t q^2 with t>0. Its translation gauge determinant is one.
Keep the original sources Jq q+Jy y. On D=1-2 v t Jy>0,

W=log Z=-log D/2+v Jq^2/(2D),
<q>=v Jq/D,
<q^2>=v/D+<q>^2,
<y>=t<q^2>.

The Euclidean Legendre convention is Gamma=J phi-W. For
r=<y>/t-<q>^2>0,
Gamma=<y>/(2vt)-1/2+log(v/r)/2.
Its Nielsen vector is K=(0,<y>/t), and
partial_t Gamma+K dot grad Gamma=0 exactly. The stationary
mean is (0,tv), while the physical q distribution is unchanged.
The full stationary Hessian is diag(1/v,1/(2t^2v^2)).
All source and Hessian identities, actual Gaussian moments and
independent numerical integrals are checked. The Legendre claim
does not extend through t=0.

A separate boundary-weight test multiplies by exp(-eta y^2/2).
After the same reduction the physical q weight contains
exp(-eta t^2 q^4/2). Its normalized <q^2> has second derivative
-12 eta v^3 at t=0, computed by Gaussian moments and independent
integration. This limit concerns the physical integral, not the
singular two-observable Legendre chart. A gauge-dependent
normalization at zero source does not remove the physical change.
