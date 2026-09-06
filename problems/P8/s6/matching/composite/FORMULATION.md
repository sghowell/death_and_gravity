# S6.6.COMPOSITE: a different physical metric and matter hypothesis

This gate constructs a regular composite-metric free-scalar bounce and an
exact local CD history with a reconstructed canonical potential. It is
not a parent-matching or EFT-health certificate. It is a new sibling of
the unchanged HR screens, not an extension of their exclusions to shared
matter or another physical frame.

## Prescribed action and physical variables

Use the P8 `+---` convention and the same HR gravitational action as the
pinned S6.5 lineage, with `G=M_g^2>0`, `F=M_f^2>0`, `m4>0` and arbitrary
constant real `beta0,...,beta4`. Replace all old separate matter sectors
by exactly one canonical scalar:

\[
S_\chi=\int\sqrt{|g_e|}\left[\tfrac12g_e^{\mu\nu}
\partial_\mu\chi\partial_\nu\chi-V(\chi)\right],\qquad
g_e=g(\alpha I+\beta\sqrt{g^{-1}f})^2,
\quad\alpha,\beta>0.
\]

The prescribed physical metric is `g_e`, with its proper time `T`. In
common spatially flat positive-root FLRW,
`a_e=alpha*a+beta*b`, `N_e=alpha*N_g+beta*N_f`. Both scale factors and
lapses are smooth, positive and finite on the interval considered.
Only the composite source is separately conserved; the induced `g` and
`f` sources are not separately conserved. Positive kinetic energy implies
`rho+p=(D_e chi)^2>=0`. Nonnegative potential is not assumed.

## Exact claims

With `y=b/a`, `c=N_f/N_g`, `r=alpha+beta*y`, `s=alpha+beta*c`, set

\[
P=m4(\beta_1+2\beta_2y+\beta_3y^2),\qquad
Q=P-\alpha\beta r^2p,\qquad B=N_g\dot b-N_f\dot a.
\]

The full source-aware equations imply `QB=0`, without dividing either
factor. On an open dynamical branch `B=0`, the normalized Hubble
`Z=H_g/sqrt(G+F*y^2)` is nonincreasing and `H_e=H_g/r`; contraction
cannot turn into expansion within that branch. A nondegenerate bounce
at `B=0` is excluded even at `Q=B=0`. Thus a CD bounce must instead
have `B!=0` and a pressure-branch neighborhood `Q=0`.

The pressure branch admits an explicit regular bounce with `V=0` and
nonzero scalar kinetic energy. Its pressure equation has a double root
at the bounce, but the time-dependent metric solution is regular; an
inverse `y(p)` chart is not used. This free example is not the CD history.

For `G=F=M^2`, `m4=M^2/tau^2`, `alpha=beta=1`, and
`beta_n=(0,0,1,0,0)`, the regular reconstruction ODE in the
[proof](notes/background.md) supplies an exact solution with

\[
a_e(T)=[1+(T/\tau)^2]^2,\qquad |T|\le\tau/64.
\]

The proof constructs one analytic potential on the field interval
traversed by the solution. It keeps
`3/4<y<5/4`, `1/4<rho_bar<3/4`, `1/7<N_g,N_f<6/7`, and
`rho_bar+p_bar>1/2`. This is not a solution for a potential supplied in
advance, nor a claim about full-CD free M1 matter.

## Evidence and explicit limits

The [certificate](certificates/composite-background.json) replays the
unchanged prior lineage, exact source variations, both Bianchi factors,
branch intersections, free-bounce and reconstruction identities, and a
separate coefficientwise Fraction calculation. The analytic ODE and
inverse-function arguments are written proofs supported by exact rational
margins; they are not a numerical-integration or formal PDE certificate.
The [source audit](notes/sources.md) separates background algebra from
the generic BD-mode and cutoff issues of the bare composite coupling.

No global CD continuation, unique global potential, old CD/M1 action
matching, rolling perturbation stability, interaction/loop control,
background-dependent cutoff, positivity closure or healthy UV completion
is claimed. A known trimetric higher-source completion is a different
action; its corrections are not silently omitted or transferred here.
General S6 matching and P8 remain open.
