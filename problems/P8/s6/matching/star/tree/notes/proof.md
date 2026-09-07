# Finite HR-tree no-bounce proof with arbitrary algebraic strata

## 1. Literal edge action and normalization

Adopt the action, regular actual metrics and separately conserved isotropic
NEC matter in [FORMULATION.md](../FORMULATION.md). For one oriented edge
\(i\to j\), put \(y=a_j/a_i>0\), \(c=N_j/N_i>0\), and

\[
 U=\beta_0+3\beta_1y+3\beta_2y^2+\beta_3y^3,\quad
 V=\beta_1+3\beta_2y+3\beta_3y^2+\beta_4y^3,
\]

\[
 P(y)=2(\beta_1+2\beta_2y+\beta_3y^2).
\]

The lapse-retaining edge density is \(L_e=-2N_i a_i^3(U+cV)\).
Equivalently its elementary symmetric polynomial generating function is
\((1+zc)(1+zy)^3\). This fixes every multiplicity and the factor two.
S6.5's convention is recovered by replacing its \(m^4\beta_n\) with
the present \(2\beta_n\); no numerical coefficient is imported without
that dictionary. P8(b) uses \(+---\), \(R_B=-6(DH+2H^2)\), and
EH \(-G R_B/2\), not the opposite A/FK curvature convention.

Literal lapse and scale variations give
\(\rho_{ie}=-a_i^{-3}\partial L_e/\partial N_i\) and
\(\pi_{ie}=(3N_i a_i^2)^{-1}\partial L_e/\partial a_i\), and likewise
at \(j\). They imply

\[
 \rho_{ie}=2U,\qquad \rho_{je}=2V/y^3,
\]

\[
 I_{ie}:=\rho_{ie}+\pi_{ie}=(y-c)P,
 \qquad I_{je}:=\rho_{je}+\pi_{je}={c-y\over cy^3}P.
 \tag{1}
\]

The null-stress reciprocity is therefore

\[
 N_i a_i^3 I_{ie}+N_j a_j^3 I_{je}=0.
 \tag{2}
\]

Endpoint coefficients contribute cosmological densities but cancel in
these null stresses. Edge reversal swaps \(\beta_n\leftrightarrow
\beta_{4-n}\) and leaves the actual action unchanged.

## 2. Source-aware Bianchi flux and finite-tree elimination

Define \(D_i=N_i^{-1}d/dt\), \(H_i=D_i\log a_i\), and the interaction
balance at one endpoint
\(C_{ie}=D_i\rho_{ie}+3H_iI_{ie}\). Before imposing any branch,
\(D_iy=y(cH_j-H_i)\), so direct differentiation of the densities gives

\[
 C_{ie}=3Pc(yH_j-H_i),\qquad
 C_{je}=-{C_{ie}\over c^2y^3}.
 \tag{3}
\]

The reciprocal **Bianchi** weight contains two lapses:

\[
 \mathcal F_e:=N_i^2a_i^3 C_{ie}
 =3N_iN_j a_i^2P(a_jH_j-a_iH_i)
 =-N_j^2a_j^3 C_{je}.
 \tag{4}
\]

This is distinct from the one-lapse null weight (2). Using it in place
of the two-lapse weight leaves the generically nonzero defect
\(3P(c-1)(yH_j-H_i)\).

At every vertex let \(\rho_v,\pi_v\) now denote total matter plus
incident interaction density and pressure. The lapse/spatial residuals are

\[
 E_{0v}=3G_vH_v^2-\rho_v,\qquad
 E_{sv}=G_v(2D_vH_v+3H_v^2)+\pi_v.
\]

The exact Bianchi identity is

\[
 D_v E_{0v}+3H_v(E_{0v}-E_{sv})
 +D_v\rho_v+3H_v(\rho_v+\pi_v)=0.
 \tag{5}
\]

It remains true when \(G_v=0\): the full vertex equations are then
algebraic, but are not discarded. On actual solutions, separately
conserved matter at that vertex gives \(\sum_{e\ni v}C_{ve}=0\).
Multiplying by \(N_v^2a_v^3\) turns these equations into

\[
 B_{\mathcal T}\mathcal F=0,
 \tag{6}
\]

where the oriented incidence matrix is \(+1\) at an edge's source and
\(-1\) at its target. A leaf has one incident flux, so that flux must
vanish. Reciprocity removes the same edge at its neighbor, and repeated
leaf stripping exhausts a finite tree. Thus every \(\mathcal F_e=0\).
Equivalently its incidence matrix has trivial column kernel. All weights
in (4) are finite and nonzero, giving the undivided edge relation

\[
 P_e(y_e)(a_jH_j-a_iH_i)=0.
 \tag{7}
\]

No Hubble rate or interaction polynomial has been divided out at a zero.
Independent rational Gaussian elimination checks the incidence result for
chains, stars and a branched tree using a different algorithm. An oriented
triangle has a nonzero circulating flux in the incidence kernel, so the
tree proof is not valid for cycles. This is a graph countercontrol, not a
full cyclic background or a proof of a cyclic bounce.

## 3. The instantaneous active component

Fix any time \(T_*\) and physical vertex \(r\), with \(G_r>0\).
Let \(C_*\) be the connected component of \(r\) after keeping precisely
the edges with \(P_e(y_e(T_*))\ne0\). There are finitely many such
internal edges. By continuity each remains nonzero in a neighborhood of
\(T_*\); their intersection is a neighborhood. Equation (7) therefore
implies there, along every internal path,

\[
 H_i=H_r/y_i,\qquad
 y_i=a_i/a_r>0,\quad c_i=N_i/N_r>0,
 \quad D_r y_i=(c_i-y_i)H_r\quad(i\in C_*).
 \tag{8}
\]

These \(c_i\) are relative lapses, not physical cone speeds. The physical
Einstein cone speed would be \(c_i/y_i\); no cone condition is used.
Every edge leaving \(C_*\) has \(P_e(T_*)=0\), so both null stresses
in (1) vanish **at that point**. It need not stay algebraic on a
neighborhood; no derivative of its zero condition is used.

This observation covers persistent algebraic roots, isolated roots,
multiple roots, arbitrarily accumulated switches and identically zero
polynomials without a root-set classification. For a genuine polynomial,
a persistent algebraic ratio is indeed constant on any connected root
interval, but that fact is unnecessary for the present proof. Endpoint-only
edges simply remain inactive. The root always lies in its own component,
even if that component is a singleton.

## 4. Conditional null combination with positive target inertia

The vertex null equation is
\(-2G_iD_iH_i=n_i+\sum_e I_{ie}\), with actual matter
\(n_i\ge0\). At \(T_*\), multiply it by

\[
 w_i=c_i y_i^3={N_i a_i^3\over N_r a_r^3}>0
\]

and sum over \(i\in C_*\). Internal terms cancel by (2); boundary
terms vanish by their pointwise \(P_e=0\). Substituting (8) and
\(D_i=c_i^{-1}D_r\), valid for the fixed internal component near
\(T_*\), yields

\[
 -2A_C D_rH_r+H_r(D_rA_C)_{\rm fixed\ C}
 =\mathcal N_C,\qquad
 A_C=\sum_{i\in C_*}G_i y_i^2\ge G_r>0,
 \quad \mathcal N_C=\sum_{i\in C_*}c_i y_i^3n_i\ge0.
 \tag{9}
\]

Only the vertex subset selected at \(T_*\) is held fixed while
differentiating. The result does not differentiate a component label that
can jump as time varies. It does not assert a global exact identity for
\(H_r/\sqrt{\sum_{\rm all}G_i y_i^2}\).

Define at the point

\[
 \lambda_C={\sum_{i\in C_*}G_i y_i^2 D_r\log y_i\over A_C},
 \qquad \eta_C={\mathcal N_C\over2A_C}\ge0.
\]

Then \(D_rH_r=\lambda_C H_r-\eta_C\). Since \(G_i\ge0\),

\[
 |\lambda_C|\le\max_{i\in C_*}|D_r\log y_i|
 \le M(T):=\max_{i\in\mathcal T}|D_r\log(a_i/a_r)|.
 \tag{10}
\]

There are finitely many positive smooth ratios. On every compact regular
physical-time interval, \(M\) is continuous and bounded. The strictly
positive target coefficient supplies a denominator bound regardless of
which edges become algebraic or which other Einstein coefficients vanish.
At every point, including all algebraic boundaries,

\[
 D_rH_r\le M(T)|H_r|.
 \tag{11}
\]

## 5. All-point positive-part proof and strict signs

Choose an arbitrary compact ordered subinterval \([T_0,T_1]\) in the
regular domain, and let \(M_*=\sup M<\infty\). The positive part
\(Y=\max(H_r,0)\) is absolutely continuous. Where \(H_r>0\), (11)
gives \(Y'\le M_*Y\). Where \(H_r<0\), \(Y'=0\). On the zero set
the derivative of the positive part is zero almost everywhere; moreover
(9) gives \(H_r'\le0\) at every zero. Hence

\[
 Y'\le M_*Y\quad\text{almost everywhere},\qquad
 (e^{-M_*(T-T_0)}Y)'\le0.
\]

If \(H_r(T_0)\le0\), then \(Y(T_0)=0\), and nonnegativity forces
\(Y=0\) throughout. Exhausting compact subintervals proves the stated
global sign result on the connected regular interval, with no assumption
on a tail limit or a finite number of branch switches.

If \(H_r(T_0)<0\), then up to any putative first zero, \(W=-H_r>0\)
satisfies \(W'\ge-M_*W\), hence
\(W(T)\ge W(T_0)e^{-M_*(T-T_0)}>0\). Continuity rules out that zero.
Thus strict negativity also persists. These arguments exclude every
degenerate negative-to-positive crossing, not only a point with strictly
positive acceleration. They do not forbid expansion or a turnaround.

An optional exact integrating-factor formulation sets
\(\lambda(T)=\lambda_{C(T)}(T)\). The finite-valued component selection
is Borel measurable, and (10) makes \(\lambda\) locally bounded. Then
\(e^{-\int\lambda dT}H_r\) is absolutely continuous and nonincreasing,
with derivative \(-e^{-\int\lambda dT}\eta_{C(T)}\) almost everywhere.
This is a diagnostic integrating factor, not a local physical field map,
a canonical normalization, or the derivative of a changing \(A_C\).
The simpler proof above does not require using it.

## 6. Actual examples, omissions and target-kinetic guard

The actual expanding controls have nonzero separately conserved sources,
not merely prescribed instantaneous null densities.

For any finite tree with aligned metrics
\(a_v=T^{1/3},N_v=1,T>0\), put
\(\beta_e=(-3,1,0,0,-1)\) on every edge. All interaction stresses
vanish at the common ratio/lapse one, while \(P_e=2\ne0\). Give each
vertex its own free scalar trajectory

\[
 \phi_v=\sqrt{2G_v/3}\log T,\qquad V_v=0.
\]

For \(G_v=0\) the field is constant. At each vertex
\(H_v=1/(3T)\), \(\rho_v=p_v=G_v/(3T^2)\), both Einstein equations
and the conserved free-scalar current hold. The checked four-vertex
branched example uses \(G=(1,0,2,3)\), including a zero-Einstein internal
vertex, and three orientations. These are local regular expanding
solutions, not bounces, spectrum or cutoff witnesses.

A three-vertex mixed algebraic control has \(G_v=1\), equal
\(a_v=T^{1/3}\), lapses \((N_0,N_1,N_2)=(1,1/(3T),1)\), and
oriented edges

\[
 1\to0:\ \beta=(0,1,-1/2,0,1/2),\qquad
 2\to0:\ \beta=(-3,1,0,0,-1).
\]

Vertices 0 and 2 have independent free scalars
\(\sqrt{2/3}\log T\), and vertex 1 has a constant scalar. The first
edge has \(P=0\); its leaf is de Sitter with \(H_1=1\), and its central
interaction density/pressure vanish. The second edge has \(P=2\).
All lapse, scale and scalar equations hold. For physical vertex 0,
\(C=\{0,2\}\), \(A_C=2\),
\(\mathcal N_C=4/(3T^2)\). Equation (9) holds, whereas replacing
\(A_C\) by the all-vertex value 3 produces residual \(2/(3T^2)\).
This is an actual action-consistent omission control, not just a sampled
coefficient choice.

The positive target coefficient cannot simply be removed. With two
vertices, \(G_0=0,G_1=1\), all edge beta and sources zero,
\(a_0=1+T^2,a_1=1,N_0=N_1=1\), every field equation holds. The actual
target has \(H_0'(0)=2\) but is undetermined, and its instantaneous
component has \(A_C=0\). Adding a positive \(G_0\) changes its equations.
This does not refute the stronger auxiliary-center star theorem under its
different genuine-link conditions; it explains the present simple target
guard.

NEC is also material. Aligned metrics with \(a=1+T^2\), \(G=1\),
tuned beta1 edges and separate fluid profiles
\(\rho=3H^2,p=-2H'-3H^2\) solve all background equations and the
conservation law. At the bounce \(\rho+p=-4\). Near the bounce these
can be reconstructed by negative-kinetic-sign scalars, not by the positive
canonical matter of the theorem. Cyclic incidence circulations and
nonzero-divergence flux solutions are only algebraic omission controls;
they are not claimed as actual cyclic or shared-matter bounces.

## 7. Conditional CD matching and research boundary

Suppose a proposed construction supplies an actual solution of this exact
tree action and identifies the actual metric/proper time of a
positive-Einstein physical vertex with the original CD target
\(H_{CD}(T)=4T/(\tau^2+T^2)\). At \(T=\pm L\tau\), errors smaller
than \(4L/[\tau(1+L^2)]\) at both endpoints force the forbidden sign
change. Therefore

\[
 \max_\pm|H_r(\pm L\tau)-H_{CD}(\pm L\tau)|
 \ge {4L\over\tau(1+L^2)}.
\]

At \(L=1/2\) the threshold is \(8/(5\tau)\); equality is not itself
excluded. A connected tuned-beta1 static Minkowski tree with zero sources
attains equality, since every edge stress and Einstein tensor vanishes.
The comparison requires the actual metric/time/endpoint dictionary, not
only an approximate light-mode equation. It neither supplies the original
DHOST/free-matter operator match nor assumes a physical vacuum-to-bounce
trajectory. Changed actions can change the equations.

The source-hashed replay and independent action/clock/graph tests support
this written exact algebra and real-analysis theorem. It is not a
proof-assistant formalization or a finite numerical sampling of arbitrary
branch histories. Primary sources already contain the reciprocal-link
and leaf-stripping method in related multimetric settings; see the
[source audit](sources.md). No general novelty, universal ghost/UV
exclusion, infinite-graph theorem, perturbative health or cutoff is claimed.
The completed scoped photon objective and frozen linear classification
are unchanged. Original P8 remains open.
