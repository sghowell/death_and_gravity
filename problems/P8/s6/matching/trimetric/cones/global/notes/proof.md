# Actual auxiliary-metric Hubble monotonicity

## 1. Literal action, clocks and full equations

Use the action and regular common spatially flat homogeneous chart of
[FORMULATION.md](../FORMULATION.md). In particular \(G_i>0\), all
coefficients are constant, at least one \(p_i\ne0\), and matter acts on
the actual \(h=u^T\eta u\) only. The P8(b) curvature convention is
\(+---\), \(R_B=-6(DH+2H^2)\), \(G_{B00}=3H^2\), with EH
\(-G_iR_B/2\). An opposite A/FK sign dictionary cannot be substituted.

For each Einstein metric, retaining its lapse before variation gives the
minisuperspace density, up to the EH time boundary,

\[
 -{3G_i a_i\dot a_i^2\over n_i}
 -2p_i(n_i a_u^3+3n_u a_u^2a_i)-2b_i n_i a_i^3.
\]

These lapse/scale equations agree with the literal covariant coframe
Euler map, not merely a gauge-fixed density:

\[
 3G_iH_i^2=2p_iR_i^3+2b_i,\qquad
 G_i(2D_iH_i+3H_i^2)={2p_iR_i^3\over c_i}+2b_i,
 \tag{1}
\]

\[
 D_iH_i={p_iR_i^3\over G_i}(c_i^{-1}-1),\qquad
 D_i={d\over d\tau_i},\quad {d\tau_i\over dT}={c_i\over R_i}>0.
 \tag{2}
\]

The endpoint coefficients \(b_i\) are not omitted. The full isotropic
auxiliary equations are, with \(P_i=p_i/R_i\),

\[
 B+3\sum_iP_i=-{\epsilon\rho\over2},\qquad
 B+\sum_iP_i(c_i+2)={\epsilon p\over2}.
\]

Their undivided null combination is

\[
 \sum_iP_i(c_i-1)={n_h\over2},\qquad n_h=\epsilon(\rho+p)\ge0.
 \tag{3}
\]

Here \(p_i\) is a link coefficient and \(p\) is the physical pressure.
The source may be general isotropic matter for (3); the later tensor and
explicit scalar controls impose additional matter assumptions.

Units are \([G_i]=[K]=\mathrm{mass}^2\),
\([p_i]=[b_i]=[B]=[n_h]=\mathrm{mass}^4\),
\([H_u]=\mathrm{mass}\), and \([K']=\mathrm{mass}^3\).
Ratios and cones are dimensionless. No Einstein coefficient is silently
identified with the old CD/M1 normalization.

## 2. Bianchi lock without dividing a Hubble rate

The effective interaction density and pressure in metric \(i\) are

\[
 \rho_i=2p_iR_i^3+2b_i,\qquad
 \pi_i=-2p_iR_i^3/c_i-2b_i.
\]

Before using an equation of motion, kinematics gives
\(D_iR_i=R_i[(R_i/c_i)H_u-H_i]\), and consequently

\[
 D_i\rho_i+3H_i(\rho_i+\pi_i)
 ={6p_iR_i^3\over c_i}(R_iH_u-H_i).
 \tag{4}
\]

To see why its left side vanishes on actual solutions, define the lapse
and spatial residuals from (1), \(E_{0i}\) and \(E_{si}\). The exact
identity

\[
 D_iE_{0i}+3H_i(E_{0i}-E_{si})
 +D_i\rho_i+3H_i(\rho_i+\pi_i)=0
\]

is the Einstein Bianchi identity in this chart. It follows from both
equations in (1), including at \(H_i=0\). We do not impose interaction
conservation as an unrelated additional source law, or divide a
differentiated Friedmann equation by \(H_i\).

For each nonzero constant link, the positive finite factors in (4) imply

\[
 H_i=R_iH_u,\qquad R_i'=R_i(1-c_i)H_u.
 \tag{5}
\]

If \(p_i=0\), (4) imposes no lock. Such a metric is omitted from the
nonempty index set \(I\). Link signs never enter the implication in (5).
No sum of links is divided out; isolated or persistent zeros of
\(\sum_iP_i\) leave this calculation intact.

## 3. Exact monotonicity and every sign-crossing degeneracy

Put \(K_i=G_i/R_i^2>0\) for \(i\in I\) and \(K=\sum_I K_i>0\).
Converting (2) by (5) to the actual physical proper time gives

\[
 E_i=K_i[H_u'+(1-c_i)H_u^2]-P_i(1-c_i)=0,
 \qquad K'=-2H_u\sum_I K_i(1-c_i).
 \tag{6}
\]

In particular the lapse conversion has not been dropped: the original
Raychaudhuri residual multiplied by \(G_ic_i/R_i^4\) is \(E_i\).
Let \(U=\sum_I P_i(c_i-1)-n_h/2\). For arbitrary values of these
residuals, an exact combination is

\[
 2\sum_I E_i-2U=2K H_u'-H_uK'+n_h.
 \tag{7}
\]

Actual equations give the advertised result

\[
 2K H_u'-H_uK'=-n_h,\qquad
 Z'=-{n_h\over2K^{3/2}}\le0,\qquad Z={H_u\over\sqrt K}.
 \tag{8}
\]

The quotient is well-defined and differentiable at every point of the
connected regular interval. Its denominator is positive from the Einstein
coefficients and finite positive ratios, not from a kinetic Hessian of
the auxiliary field. Integrating the derivative inequality between any
two finite ordered physical times proves \(Z(T_+)\le Z(T_-)\).
If \(H_u(T_-)<0<H_u(T_+)\), positivity of \(K\) would give the opposite
strict ordering. Thus no contraction-to-expansion crossing exists.

This proof includes a zero with \(H_u'=0\), any higher-order degeneracy,
an interval of zeros, and a nonanalytic smooth crossing. It does not need
a polynomial expansion at the zero. It does not forbid expansion,
contraction, or expansion-to-contraction. No uniform positive lower bound
on \(K\) at an infinite endpoint, completeness of either Einstein proper
clock, asymptotic limit, tensor cone, vacuum, or TT inverse is used.
The smoothness needed is simply that the displayed full equations and
Bianchi/chain derivatives exist; smooth coframes and source suffice.

At \(H_u=0,H_u'>0\), (7) is already strictly positive for \(n_h\ge0\).
For the kinematic diagnostic \(K=1,H_u=T^3\), it would require
\(n_h=-6T^2\); this is a NEC-violation diagnostic, not a constructed
parent solution.

## 4. Quantitative actual-background tests

Suppose an explicitly supplied physical matching map identifies this
parent's actual \(h\) and its proper \(T\) with the CD target

\[
 a_{CD}(T)=(1+(T/\tau)^2)^2,\qquad
 H_{CD}(T)={4T\over\tau^2+T^2},\quad \tau>0.
\]

At \(T=\pm L\tau\), \(L>0\), the absolute target Hubble value is
\(4L/[\tau(1+L^2)]\). If the maximum of the two actual parent errors is
strictly smaller, the parent Hubble values have forbidden negative and
positive signs. Equivalently a necessary condition for such a matched
actual solution is

\[
 \max_{\pm}|H_u(\pm L\tau)-H_{CD}(\pm L\tau)|
 \ge {4L\over\tau(1+L^2)}.
 \tag{9}
\]

At \(L=1/2\) this is \(8/(5\tau)\), and at \(L=1\) it is
\(2/\tau\). Equality does not itself contradict (8). No Hubble-derivative
error estimate is required. This is stronger than a full-parent-cone
screen: even an asymmetric light-only construction cannot claim a
controlled actual parent background with errors below (9), under the
specified metric/proper-time/endpoint dictionary. A calculation in an
approximate reduced equation alone supplies no such dictionary or parent
solution. The geometric test is not an operator or matter-content match,
and no physical vacuum-to-bounce trajectory is assumed.

A separate generic residual threshold makes omitted terms explicit. If
\(K\ge\kappa>0\), \(H_u'\ge a>0\), \(n_h\ge\nu\ge0\),
\(|H_u|\le\eta\), \(|K'|\le D\), then

\[
 2K H_u'-H_uK'+n_h\ge2\kappa a+\nu-\eta D.
 \tag{10}
\]

A positive lower bound requires at least that cancellation by a proposed
correction to this actual combined equation. It is not a computed bound
on loops, individual operators, field redefinitions, connection errors or
the correction to \(K\) itself. An action with extra operators can change
(4) as well as (6); this gate does not apply unchanged to that action.
Numeric interfaces use exact finite rationals with explicit signs;
ordinary floats, booleans, nonfinite and nonreal values are rejected.

## 5. Independent, weaker cone/three-slice proof

For an alternative check, impose the extra individual Einstein cone
condition \(c_i\le1\). This section is not needed for (8).

If some \(p_i=-s_i<0\), the lapse equation in (1) forces

\[
 b_i=s_iR_i^3+{3G_iH_i^2\over2}>0,\qquad
 R_i\le R_{\rm cap}=(b_i/s_i)^{1/3}.
\]

Thus \(b_i\le0\) admits no regular real flat-FLRW point at all. For
\(b_i>0\), (2) and \(c_i\le1\) imply \(D_iH_i\le0\). The function
\(\log a_i\) is concave in ordered Einstein proper time. Meanwhile
\(a_i\ge a_u/R_{\rm cap}\). If \(a_u\) grows unboundedly in both
physical time tails, choose two finite endpoints whose lower bounds on
\(a_i\) exceed its value at a fixed middle slice. For positive left and
right proper durations \(\ell,r\), concavity gives

\[
 {\log a_i(0)-\log a_i(-)\over\ell}
 \ge {\log a_i(+)-\log a_i(0)\over r},
\]

whereas the chosen endpoints make the left side negative and the right
side positive. Equal proper durations or infinite parent proper time are
unnecessary. Neither a cap nor this proof requires an actual flat vacuum.

A finite CD three-slice sufficient condition is
\(r_0(1+L^2)^2>1\), where
\(r_0=R_i(0)/R_{\rm cap}\in(0,1]\). If relative physical scale errors at
all three slices are at most \(\epsilon_a<1\), replace it by

\[
 r_0(1+L^2)^2(1-\epsilon_a)>1+\epsilon_a.
\]

The resulting strict error threshold is
\([r_0(1+L^2)^2-1]/[r_0(1+L^2)^2+1]\) when positive. A failed or equality
test is inconclusive. Here \(r_0\) is a cap ratio; it is a vacuum ratio
only after the separate calibration below. This scale test differs from
the stronger Hubble endpoint test (9).

If instead all links are nonnegative, (3), NEC and \(c_i\le1\) force
\(n_h=0\) and \(c_i=1\) for every nonzero link. Equation (2) makes each
linked \(H_i\) constant, and its constraint fixes its positive \(R_i\)
constant. Then \(H_u=H_i/R_i\) is constant and cannot produce two growing
tails. Disconnected fields do not change this implication.

## 6. Same-action vacuum and the separate tensor inverse corollary

Assume additionally an actual constant positive Lorentz-flat vacuum of
the same action, with zero scalar gradient and appropriate constant
potential. The full linked coframe equation, not just its lapse,
forces proportionality \(e_{i0}=u_0/R_{i0}\), \(R_{i0}>0\), and

\[
 b_i=-p_iR_{i0}^3,\qquad
 B_{\rm eff}=-3S_0,\quad
 B_{\rm eff}=B+\epsilon V_0/2,\quad
 P_{i0}=p_i/R_{i0},\quad S_0=P_{g0}+P_{f0}.
\]

The scalar equation must also hold, e.g. \(V'(\psi_0)=0\). Satisfying the
Einstein calibration alone does not prove this full vacuum exists.

For a nonzero \(S_0\), use the actual vacuum \(u_0\) volume/clock and
coframe perturbations \(e=(I+H/2)/R_{g0}\),
\(v=(I+J/2)/R_{f0}\). Literal stationary elimination gives the density

\[
 -2S_0\det\!\left(I+{P_{g0}H+P_{f0}J\over2S_0}\right)
 +2P_{g0}\det(I+H/2)+2P_{f0}\det(I+J/2).
\]

Its full traceful quadratic part is

\[
 {P_{g0}P_{f0}\over4S_0}
 \{[\operatorname{tr}(H-J)]^2-\operatorname{tr}[(H-J)^2]\}.
 \tag{11}
\]

This source-action bridge precedes the coefficientwise determinant test;
it is not a FP claim inferred from TT alone. Independently authored
noncommuting traceful coframe fixtures check the same invariant. Physical
vacuum kinetic coefficients are \(K_{i0}=G_i/R_{i0}^2>0\), so

\[
 q_{\rm eff,0}={2P_{g0}P_{f0}\over S_0},\qquad
 m_{FP,0}^2=q_{\rm eff,0}(K_{g0}^{-1}+K_{f0}^{-1}).
\]

For genuine mixed signs, finite positive \(m_{FP,0}^2\) requires
\(S_0<0\). At any regular same-action flat-FLRW point, (1) now implies
the sign-free calibrated identity

\[
 P_i-P_{i0}
 =-{3G_iH_i^2\over
 2R_iR_{i0}(R_i^2+R_iR_{i0}+R_{i0}^2)}\le0.
 \tag{12}
\]

Thus \(S=\sum_iP_i\le S_0<0\): the algebraic TT denominator cannot
vanish or cross zero anywhere on this specified mixed-sign branch. The
code checks (12) as a genuine residual identity with (1), not by adding a
formula to its own negative.

For homogeneous canonical scalars with no independent TT matter response,
literal physical volume/clock pullback gives
\(K_i^{TT}=G_i/(c_iR_i^2)>0\),
\(F_i^{TT}=G_ic_i/R_i^2\). Nonsingular auxiliary TT elimination adds only
algebraic relative potential and probe contacts; both full TT principal
squared speeds remain \(c_i^2\). Therefore in this additional calibrated
case a physical full-TT cone statement can be reduced to the individual
Einstein cones. It says nothing about scalar/vector health or a finite-band
light-only approximation. Positive algebraic vacuum mass is not an
invertible rolling Green operator. The main background theorem used
neither this vacuum nor this TT inference.

## 7. Genuine controls and nonempty regular solutions

The audit distinguishes full solutions from kinematic or algebraic points.

- A full mixed singular flat solution has \(p_g=1,p_f=-1,B=0\),
  \(b_g=-1,b_f=1\), \(e=v=u=I\), and zero source. It has \(S=0\).
  The TT auxiliary equation imposes \(\gamma_g=\gamma_f\), not an
  inverse or a finite FP formula. It is covered by (8), and lies outside
  the inverse corollary of section 6.
- If all links and endpoints vanish, with \(B=V=0\) and a constant free
  scalar, \(e=v=I\) and \(u=(1+t^2)I\) satisfy all equations. Actual
  proper time is \(T=t+t^3/3\), and \(H_u'(0)=2\). This is a genuine
  but underdetermined arbitrary physical metric, with no auxiliary tensor
  inverse. It explains the nonempty-link hypothesis, not a viable matching
  parent that evades it.
- With \(G_i=p_i=1,B=-6,b_i=-5/8\), constant free scalar and
  \(e=v=u\) a common de Sitter metric of \(H=1/2\), all equations hold.
  The theorem allows this zero-NEC expansion; it does not require a static
  spacetime or forbid all cosmology.
- For \(p_g=-2,p_f=1,b_g=2,b_f=-1,B=4\), identity coframes satisfy the
  calibrated Einstein equations but fail the \(u\) equation. Positive
  cap/calibration data do not create a full vacuum. Importing the cap from
  a different endpoint action is also rejected.

There is also a nonempty local free-scalar rolling family for both link
signs, not merely an algebraic point. In dimensionless units set
\(G_i=\epsilon=1,p_i=q,B=-6q,b_i=-q,V=0\), let \(e=v=r\), and
\(u=\operatorname{diag}(N,A,A,A)r\). In common \(r\) proper time,

\[
 N={A\over6A-5},\quad \rho=p={12q(A-1)\over A},\quad
 H_r^2={2q(A^3-1)\over3},\quad
 \dot A=-{6H_rA(A-1)\over6A-5},\quad \dot a_r=H_ra_r,
\]

\[
 \dot\psi=N\sqrt{2\rho},\qquad
 \dot H_r=-{6qA^3(A-1)\over6A-5},\quad
 H_u=H_r/A,\quad dT=Ndt.
\]

Choose the positive root for \(H_r\). The open domain
\(A>5/6\), \(q(A-1)>0\), \(a_r>0\) makes the ODE analytic, the
source positive and the clocks regular. For \(q<0\) this specifically
means \(5/6<A<1\); for \(q>0\), \(A>1\). The constraints, both
accelerations, every \(u\) component, and the conserved scalar current
\(a_h^3d\psi/dT\) all hold. Four independent exact rational fixtures
include both domains and check (8). The analytic local ODE gives actual
local solutions, not a specified uniform duration. In particular the
negative-link continuation does not assert a healthy flat massive
spectrum, global continuation, perturbative stability, or cutoff.

## 8. Evidence and closure boundary

The source-hashed report replays the immutable S6.14 certificate and its
ancestors, then checks the new exact symbolic identities, independent
Fraction/polynomial/dual-number calculations, their primary-engine bridges,
and separately authored covariant/action audits. Coefficient algebra plus
the written calculus and local-existence arguments establish the scoped
theorem; this is not a proof-assistant formalization or a finite-sampling
claim about all times.

The no-bounce statement applies to this constant beta0/beta1 auxiliary
action with optional separate endpoint terms and the stated source/frame
conditions. A different potential or derivative operator can change the
Bianchi or source equations. The proof does not import any mutable sibling
gate, require that a vacuum and bounce be joined by a physical trajectory,
or replace S6 controlled matching by an assumed full UV completion.
It supplies no loop, cutoff, anisotropic/curved, singular, non-NEC or
universal-parent verdict. This does not change the completed scoped photon
objective or frozen linear classification; original P8 remains open.
