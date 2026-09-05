# Canonical-matter causal robustness

This is a local, two-derivative **principal** theorem. No loop counterterm or
finite-frequency dispersion relation is calculated here.

## 1. Exact sharp normal form

Let `K=[[A,B],[B,Q]]` be positive definite and `G=[[C,D],[D,Q]]`, with `Q>0`.
The equality of the two matter entries is an essential frozen-M1 hypothesis.
Define

`k=A-B^2/Q>0`, `a=A-C`, `b=B-D`,

`sigma=(a-2Bb/Q)/k`, `r=b/sqrt(kQ)`.

The shear `S=[[1,0],[-B/Q,1]]` sends K to `diag(k,Q)` and G to
`[[k(1-sigma),-b],[-b,Q]]`. Positive diagonal normalization then sends K to I
and G to `[[1-sigma,-r],[-r,1]]`. Congruence preserves the generalized speeds,
strict kinetic/gradient health and the physical-cone inequality `G<=K`.
Therefore

`c^2_+,-=1-sigma/2 +/- sqrt(sigma^2/4+r^2)`.

Strict gradient health is equivalent to `1-sigma-r^2>0`: the lower-right
gradient entry is one, and its Schur complement must be positive. Health plus
subluminality is exactly

`b=0` and `0<=a<k`.

Indeed `det(K-G)=-b^2`; positive semidefiniteness forces b=0 and then a>=0,
while the gradient Schur complement gives a<k. Conversely those inequalities
give one luminal speed and another `1-a/k` in `(0,1]`.

For every nonzero b, the exact superluminal excess is

`c_plus^2-1=(sqrt(sigma^2+4r^2)-sigma)/2>0`.

At double luminality this is `abs(r)`. A fixed positive sigma reduces it to
`r^2/sigma+O(r^4)` at small r but never makes it zero. Thus diagonal buffering
cannot give an open causal neighborhood against mixing mismatch in frozen M1.
This is stronger than a claim based only on a small-matter expansion.

The characteristic-polynomial and congruence identities are exact symbolic
checks, accompanied by the preceding real-matrix lemma. They do not diagnose
a ghost from a singular finite-k Legendre chart.

## 2. Crossing-safe and covariant dictionaries

The earlier metric-derived action defines

`w=l*(3delta-1)`, `J0=Sigma_total+3Theta^2/T`, `J=J0-w^2/2`,

`E=Lambda-T*(1-3delta)`, `P_gamma=Theta^2*F_S`.

In the ordinary chart `b=lE/(2Theta)`; in the regular principal gamma chart
`b=-lE/(2Lambda)`. Both give exactly

`r^2=l^2 E^2/(2T^2 J)`,

`sigma=(T^2J0-P_gamma-TwlE)/(T^2J)`.

These invariants have neither a Theta nor a Lambda denominator. They remain
finite through the pinned bounce because T,J are positive. At any point their
physical interpretation still uses one of the proved covering regular charts.
The gamma kinetic action is the derived high-q principal action, not an
assertion of a positive finite-q Hamiltonian at every momentum.

For `A1=0`, direct substitution of the ADM coefficients gives
`E=-2X(2F2X+XA3)`. Since pinned M1 has `l!=0` at every finite time, causal
mixing requires `A3=-2F2X/X` on the background. The pinned CD functions satisfy
this identity throughout the clock tube, a stronger statement than equality
only at X=1. Tuning FXX cannot repair a nonzero E.

As a concrete negative control at the bounce, take

`K=[[6,1/20],[1/20,1/2]]`, `G=diag(599/100,1/2)`.

Both are positive definite and `K11-G11=1/100>0`, but
`c^2=1 +/- 1/sqrt(1199)`. Conversely, changing the matter diagonal can change
the conclusion: `K=I`, `G=[[9/10,1/20],[1/20,9/10]]` has speeds squared
17/20 and 19/20. This is a useful control outside frozen M1, not a permitted
escape within it.

## 3. An actual one-sided covariant deformation

Let J_star be the certified CD coefficient and choose a fixed finite epsilon.
In dimensionless units add only

`delta_F=epsilon*J_star(phi)*(X-1)^2/(4d^2)`, `d=1+phi^2`.

Its value, first X derivative and all background phi derivatives of those
zero functions vanish on X=1. The entire backreacted background, completeness,
matter equation, tensor coefficients, Theta, Lambda and delta are unchanged.
Ia degeneracy and the clock-tube coefficient denominator are unchanged because
F does not enter those relations. Direct expansion of `N*delta_F(phi,N^-2)`
gives

`delta_Sigma=2*delta_FXX=epsilon*J_star/d^2`,

`J_new=(1+epsilon/d^2)*J_star`.

On either old regular chart the kinetic update has only its 11 entry nonzero,
while G remains the pinned K. The exact characteristic polynomial factors as

`det(G-z*K_new)=det(K_star)*(1-z)*[1-(1+epsilon/d^2)z]`.

For epsilon>=0, K_new and G are strictly positive at every finite time and
the speeds squared are `{1,1/(1+epsilon/d^2)}`. For `-1<epsilon<0`, the same
matrices are healthy but one speed is superluminal at every finite time.
Epsilon=-1 loses a kinetic direction at the bounce. This is a one-sided family
inside the exceptional relation, not a mismatch robustness margin.

The d^-2 factor is chosen to preserve the canonical tails. For the old p=4
clock normalization `psi=sqrt(8)*asinh(phi)`, write `Z=(partial psi)^2`, so
`X=dZ/8`. The additional potential, kinetic and quartic coefficients are

`epsilon*J_star/(4d^2)`, `-epsilon*J_star/(16d)`, `epsilon*J_star/256`.

The exact limit `d*J_star->4` makes all three vanish. To certify the old
quantitative derivative-tail contract as well, set `P(x)=d*J_star`,
`Z=1-x^2=1/d`, `D_psi=Z*partial_x/sqrt(8)`. Each correction is
`kappa*epsilon*P*Z^n`, with `(kappa,n)=(1/4,3),(-1/16,2),(1/256,1)`.
Define polynomials `P1=Z*P'-2n*x*P`, `P2=Z*P1'-2n*x*P1`. Their first and
second canonical derivatives are exactly `kappa*epsilon*Z^n*P1/sqrt(8)` and
`kappa*epsilon*Z^n*P2/8`. The verifier records every polynomial and explicit
rational l1 majorants: the first-derivative square decays as d^(-2n) and the
second derivative as d^(-n). The three coefficient value bounds are
`2*abs(epsilon)/d^3`, `abs(epsilon)/(2d^2)` and `abs(epsilon)/(32d)` using
the independently certified P<8. X and mixed derivatives of these added
coefficient functions vanish identically; the other covariant tail terms are
unchanged.

The old canonical kinetic error is bounded by `(10481/6720)/d`. Hence, for
epsilon>=0, `d>=10481/1680+2*epsilon` is a sufficient explicit condition for
the new canonical kinetic coefficient to be at least 1/4. The finite-time
principal health is already proved by the factorization above. No uniform
tail-entry threshold independent of unbounded epsilon is claimed.

## 4. Source and scientific boundary

[Mironov, Rubakov and Volkova, arXiv:2011.14912v1](https://arxiv.org/pdf/2011.14912v1),
equations (14)--(19), supplies the canonical-matter context and exceptional
DHOST Ia subclass. Its conventions match ours, but its kinetic/gradient matrices
are named G_AB/F_AB. Its equations (17)--(18) use small matter density; the
normal form above is derived directly and exactly, with no such expansion.
The old covariant/action and crossing certificates are retained as pinned input.

The theorem does not prove radiative instability or technical naturalness.
There is no calculation showing that any loop generates a particular E or a,
or that radiative corrections stay within quadratic Ia plus precisely the
original canonical matter action. Higher derivatives, matter-curvature terms,
changed matter kinetic terms or lost degeneracy require a new reduction.
Even a calculated speed excess would need to exceed EFT truncation uncertainty
before it could be identified as a resolved effect. Perturbative control and
the separate UV assumptions remain open.
