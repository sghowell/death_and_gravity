# Proof of the specified beta1 parent screen

This is an exact test of the [new parent](../FORMULATION.md), in fixed
physical `g` variables. It has two distinct results: a positive constant
vacuum and flat quadratic source match, and an obstruction to the required
rolling CD geometry. Passing the first does not bypass the second.

## 1. Lapse-retaining reduction

Write `G=M_g^2>0`, `F=M_f^2>0`, `nu>0`, and

\[
g=N_g^2dt^2-a^2d\mathbf x^2,\quad
f=N_f^2dt^2-b^2d\mathbf x^2,\quad
y=b/a>0,\quad c=N_f/N_g>0.
\]

The eigenvalues of the specified square root are `(c,y,y,y)`. In the P8
curvature convention,

\[
R(a,N)=-6\left(\frac{\ddot a}{aN^2}
 +\frac{\dot a^2}{a^2N^2}-\frac{\dot a\dot N}{aN^3}\right).
\]

Subtracting the exact Einstein time boundary
`3G a^2 adot/N_g`, and similarly for `f`, gives the coordinate-volume
Lagrangian

\[
L=-\frac{3Ga\dot a^2}{N_g}-\frac{3Fb\dot b^2}{N_f}
 -\nu[-3N_ga^3+N_fa^3+3N_ga^2b-N_fb^3]
 +\frac{a^3K}{2N_g}-N_ga^3V.
\]

Here `K=G_AB(phi) phi_dot^A phi_dot^B` is the coordinate kinetic quadratic
form, not a lapse-dependent quantity; for the explicit two-field parent
it is `phi_dot^2+chi_dot^2`. Thus

\[
\rho=K/(2N_g^2)+V,\quad p=K/(2N_g^2)-V,\quad
\rho+p=K/N_g^2\geq0.
\]

No positivity of `V` is needed for the obstruction. In particular kinetic
mixing cannot change the sign when its full field metric is positive.
An independently tested LDL square decomposition checks this distinction.

Put `H_g=adot/(N_g a)`, `H_f=bdot/(N_f b)`, and `D_i=N_i^{-1}d/dt`.
Literal lapse-first Euler variation yields

\[
\begin{split}
3GH_g^2&=\rho+3\nu(y-1),\\
3FH_f^2&=\nu(y^{-3}-1),\\
G(2D_gH_g+3H_g^2)&=-p-\nu(3-2y-c),\\
F(2D_fH_f+3H_f^2)&=-\nu(1-(cy^2)^{-1}).
\end{split}\tag{1}
\]

Both lapse equations and both scale equations are required. Gauge fixing
the lapses before variation would discard constraints needed below.
The direct four-dimensional covariant audit checks the Einstein reduction
and all four normalized Euler expressions independently.

The effective interaction density/pressure pairs are

\[
(\rho_{gv},p_{gv})=(3\nu(y-1),\nu(3-2y-c)),\qquad
(\rho_{fv},p_{fv})=(\nu(y^{-3}-1),\nu(1-(cy^2)^{-1})).
\]

The `g` Noether identity, together with the canonical matter equation, gives

\[
D_g\rho_{gv}+3H_g(\rho_{gv}+p_{gv})
=\frac{3\nu}{N_g^2a}(N_g\dot b-N_f\dot a)=0.\tag{2}
\]

Since `nu>0`, the undivided dynamical branch is
`N_g bdot=N_f adot`. The alternative interaction polynomial is the constant
`beta1>0`, so there is no algebraic branch in this family. Equation (2),
unlike the often-used `N_f=bdot/adot` chart, remains meaningful at a bounce.
In particular `H_f=H_g/y` without dividing by either Hubble parameter.

## 2. Exact CD obstruction and robust mismatch

Choose physical `g` cosmic time only after (1). At any putative regular
stationary physical slice `H_g=0`, (2) gives `H_f=0`. The `f` lapse equation
then implies

\[
0=\nu(y^{-3}-1)=-\frac{\nu(y-1)(y^2+y+1)}{y^3}.
\]

The denominator is positive and
`y^2+y+1=(y+1/2)^2+3/4>0`, hence `y=1`. Smoothness of the branch gives
`D_fH_f=Hdot_g/c` there. Taking the null combinations of (1), multiplying
the `f` equation by its positive lapse ratio, gives

\[
-2G\dot H_g=\rho+p+\nu(1-c),\qquad
-2F\dot H_g=\nu(c-1).
\]

The interaction terms cancel exactly:

\[
\boxed{-2(G+F)\dot H_g=\rho+p\geq0.}\tag{3}
\]

No regular nondegenerate physical bounce exists in the stipulated family.
For CD, `a=(1+(t/tau)^2)^2` gives `Hdot_CD(0)=4/tau^2>0`, requiring the
forbidden null stress `-8(G+F)/tau^2`. At a parent stationary slice compared
to that target slice, the Hdot error is at least `4/tau^2`; increasing
`M*tau` cannot suppress this error relative to the background scale.

A shifted stationary slice cannot produce arbitrarily good matching on a
whole bounce window either. On `I=[-tau/2,tau/2]`, CD has endpoint values
`H_CD=±8/(5tau)` and

\[
\dot H_{CD}\geq\frac{48}{25\tau^2},\qquad
\tau^2\dot H_{CD}-\frac{48}{25}
=\frac{(1-4w)(52+12w)}{25(1+w)^2}\geq0,
\quad w=(t/\tau)^2\in[0,1/4].
\]

If a regular parent had endpoint `H` errors strictly below `8/(5tau)` and
uniform `Hdot` error strictly below `48/(25tau^2)`, its endpoint Hubble
values would have opposite signs, while its Hubble derivative would be
positive throughout `I`. The intermediate-value theorem would then give
a regular zero with positive derivative, contradicting (3). This is a
physical-cosmic-time `H,Hdot` mismatch, not a norm comparison in a changed
metric or an arbitrary time coordinate.

The algebraic control `G=F=nu=a=b=N_g=1`, `N_f=4/5`, `adot=bdot=0`,
`addot=1/10`, `bddot=2/25`, `K=-2/5`, `V=1/5` satisfies all four equations
at one instant and the differentiated Bianchi relation. It is deliberately
outside positive-metric matter, since `rho+p=-2/5`. It is not claimed to be
a solution or healthy witness; it shows why the matter-sign assumption
cannot be silently removed. Dropping either `f` equation also triggers a
nonzero omission control. No small-curvature approximation enters (3).

## 3. Complete quadratic constant-vacuum spectrum

At `f=g=eta`, the potential and both tadpoles vanish. For the relative
matrix perturbation `H=eta^{-1}(h_f-h_g)`, direct determinant/square-root
trace expansion gives

\[
\sqrt{|g|}[-3+e_1(\sqrt{g^{-1}f})-e_4(\sqrt{g^{-1}f})]
=\frac18[\operatorname{tr}(H^2)-(\operatorname{tr}H)^2]+O(h^3).
\]

This is the full Fierz--Pauli relative mass structure, not just a fit on TT
fields. Set

\[
M^2=G+F,\quad \mu^2=\nu(1/G+1/F),\quad
h_0=\frac{Gh_g+Fh_f}{M^2},\quad h_m=h_f-h_g.
\]

The inverse map is `h_g=h_0-F h_m/M^2`, `h_f=h_0+G h_m/M^2`. Their Einstein
kinetic coefficients are `M^2>0` and `GF/M^2>0`; the relative FP mass is
`mu^2>0`. The explicit canonical vacuum matter has a positive massive clock
and positive free massless M1 scalar. At constant vacuum values the matter
stress has no linear perturbation, hence no quadratic metric/matter mixing.

For one real unit-TT-norm mode, define `D=partial_t^2+k^2` and normalize the
external source by

\[
S_T=-\frac18\left[Gh_gDh_g+Fh_fDh_f+\nu(h_g-h_f)^2\right]
     +\frac14j h_g.\tag{4}
\]

For a covariant source convention `+h:T/2`, `j=2T`. A common spatial-metric
perturbation sign changes both source and field, not the following squared
coefficient. The source operator is normalized so `K h=(j,0)`.

The physical-source response has pole weights `1/M^2>0` and `F/(GM^2)>0`.
The conserved-source calculations below, not TT alone, verify the two
massless and five massive physical polarizations. This is a quadratic
Minkowski result. It is not a theorem about a rolling branch, interaction
strength, all nonlinear modes, or UV completion.

## 4. True massive integration and physical metric Schur complement

The mass-basis version of (4) is

\[
S_T=-\frac18\left[M^2h_0Dh_0+\frac{GF}{M^2}h_m(D+\mu^2)h_m\right]
     +\frac14j\left(h_0-\frac{F}{M^2}h_m\right).
\]

Integrating the actual massive mode through its own sourced equation gives

\[
h_m=-\frac{j}{G(D+\mu^2)},\qquad
\Delta S_{\rm source}=\frac{j^2F}{8GM^2(D+\mu^2)}.\tag{5}
\]

The nonlocal source-source term is part of the result; omitting it removes
physical massive exchange. The remaining physical observable is still
`h_g=h_0-F h_m/M^2`, not just `h_0`. Use a common chosen vacuum Green
prescription (for example Feynman) in all inverses. The displayed identities
are meromorphic kernel identities away from poles, not a prescription for
arbitrary homogeneous massive initial data.

Independently, eliminating `f` through its own equation while keeping `g`
as the sourced variable gives the exact spin-2 Schur kernel and response

\[
\begin{split}
K_2(D)&=GD+\nu-\frac{\nu^2}{FD+\nu}
       =M^2D-\frac{F^2D^2}{\nu+FD},\\
K_2^{-1}(D)&=\frac1{M^2D}+\frac{F}{GM^2}\frac1{D+\mu^2}.
\end{split}\tag{6}
\]

Equation (5), including the physical observable map, reproduces (6).
Nevertheless `f` is not the pure massive eigenfield: its inverse auxiliary
kernel pole at `D=-nu/F` is not the physical massive pole `D=-mu^2`.

Solving the *g* equation for `f` and substituting it in the source-free
action is a different operation. It gives a `D^2` kernel coefficient
`(G^2+2GF)/nu` instead of `-F^2/nu`, a mismatch `(G+F)^2/nu`. This is an
exact negative control against silently relabeling that action as (5).

## 5. Scalar/source projectors and the curvature-squared coefficient

It is essential to check the trace-sensitive part of exchange. At a non-null
timelike momentum a conserved source has only a spatial symmetric tensor.
Use six orthonormal symmetric components (the off-diagonals multiplied by
`sqrt(2)`). Let `P0` project onto its spatial trace, and `P2=I-P0` onto its
five traceless components. A direct all-ten-component FP variation gives
these projectors rather than assuming them from TT.

In detail, for spatial perturbation `s`, lapse perturbation `n`, mixed
components `v_i`, kinetic coefficient `kappa>0`, and FP mass `mu^2>0`,

\[
S^{(2)}=-\frac{\kappa D}{8}[s:s-(\operatorname{tr}s)^2]
 -\frac{\kappa\mu^2}{8}[s:s-(\operatorname{tr}s)^2
                       +2n\operatorname{tr}s-2v_iv_i]
 +\frac14j:s.
\]

The lapse equation sets `tr s=0`, mixed equations set `v_i=0`, and
`n=tr j/(3 kappa mu^2)`. Thus
`s=P2 j/[kappa(D+mu^2)]`. Removing the mass term and gauge fixing the
temporal components instead gives `s=(P2-P0/2)j/(kappa D)`. Consequently
the full conserved physical-source response is

\[
\mathcal G(D)=\frac{P_2-P_0/2}{M^2D}
       +\frac{F}{GM^2}\frac{P_2}{D+\mu^2},\qquad
\mathcal K(D)=K_2(D)P_2-2M^2DP_0.\tag{7}
\]

The contractions are `T:T-(tr T)^2/2` for the massless term and
`T:T-(tr T)^2/3` for the massive one. The negative off-shell `P0` coefficient
in (7) is a constraint/gauge contribution, not a negative physical residue.
On the null pole, solve conservation with momentum `(1,0,0,1)` directly:
the massless contraction reduces to `(T11-T22)^2/2+2 T12^2`. In the massive
rest frame the contraction is

\[
\frac{(T_{11}-T_{22})^2+(T_{22}-T_{33})^2+(T_{33}-T_{11})^2}{3}
 +2(T_{12}^2+T_{13}^2+T_{23}^2).
\]

These are respectively two and five positive physical squares. No
division by null momentum squared is used in the null-pole check.

The direct four-index curvature contraction for arbitrary transverse metric
components, or equivalently the flat quadratic covariant identities modulo
the Euler boundary, gives

\[
C_{\rm lin}^2=\tfrac12D^2 h:P_2:h,\qquad
R_{\rm lin}^2=3D^2 h:P_0:h.
\]

The direct calculation retains all six spatial entries at timelike
momentum. Lorentz covariance and polynomial continuation from that open
non-null domain give the flat covariant quadratic identity; this is not
an analyticity assumption about cosmological amplitudes. Since (7) expands
to

\[
\mathcal K_4=M^2D(P_2-2P_0)-\frac{F^2}{\nu}D^2P_2,
\]

the action `-h:K:h/8` is reproduced by

\[
\boxed{c_C=\frac{F^2}{4\nu}
 =\frac{M^2\alpha^2}{4\mu^2}>0,\quad
c_R=0,\quad \alpha^2=F/G.}\tag{8}
\]

Here `c_C,c_R` multiply `sqrt|g| C^2,sqrt|g| R^2` in the **flat quadratic
tree-level pure-metric representative** with prescribed source and modulo
total derivatives. Only these two four-derivative quadratic coefficients
are determined. Higher interactions, nonlinear covariant matching, a
rolling-background matching expansion, and loops are not thereby fixed.
Integrating in the mass basis leaves the source functional (5); presenting
only a massless Einstein action and discarding its source operators is not
an alternative derivation of zero physical matching effects.

## 6. Exact flat remainder and scope of scales

The exact rational remainder in (6) is

\[
K_2-K_{2,4}=\frac{F^3D^3}{\nu(\nu+FD)}.
\]

For `0<=eta<1` and `|FD|<=eta nu`, the reverse triangle inequality gives

\[
|K_2-K_{2,4}|\leq\frac{F^3|D|^3}{\nu^2(1-\eta)},\qquad
\frac{|K_2-K_{2,4}|}{M^2|D|}
 \leq\frac{\eta^2}{1-\eta}.
\]

The second estimate uses the exact identity
`(K2-K2,4)/(M^2D)=(F/M^2)r^2/(1+r)`, `r=FD/nu`, and `0<F/M^2<1`.
Its value at `D=0` is understood continuously. These are spectral-kernel
bounds; `D=-omega^2+k^2` is not `k^2`. Near a massless shell it can be
small even for large spatial momentum. No spatial-band cutoff, interaction
bound, or time-dependent Green-function error follows from this alone.
Similarly a parametrically large `mu` does not prove weak coupling.

In natural units `G,F,M^2` have mass dimension two, `nu` dimension four,
`D,mu^2` dimension two, `c_C,c_R` dimension zero, and
`beta=c_C/M^2=alpha^2/(4mu^2)` dimension minus two. All terms in (3) have
dimension four. None of these scale dictionaries repairs its sign
contradiction at arbitrarily large `M*tau`.

## 7. Certification boundary

The report replays both pinned inputs and exact new algebra; it hashes this
proof, the formulation/source notes, all source modules, and all tests,
including the separate covariant audit. Pytest also checks exact examples,
the full FP source constraint, lapse/source omissions, forbidden-input
controls, and certificate mutations. Polynomial equalities are exact, not
sampled interval evidence. The implication from positivity, branch
smoothness and the intermediate-value theorem is a written argument.

The outcome is rejection of the stipulated parent as a CD realization.
It does not reject general bimetric models, quantum matter, changed physical
frames, other operators, non-flat geometries, or UV completions. It does not
turn the earlier conditional local-Weyl candidate into a completed parent
match. A zero flat tree R² coefficient does not cancel the mandatory
isolated-matter one-loop R² running, unknown finite coefficients, anomaly
terms, or massless nonlocal contributions. No exact fourth-order resummation
or causal/ghost inference from a truncated curvature action is made.
