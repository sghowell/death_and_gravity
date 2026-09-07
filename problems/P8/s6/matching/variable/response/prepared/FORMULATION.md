# Prepared analytic physical tensor sector

This descendant studies the **unchanged** S6.20 VARIABLE action and actual
local background, through frozen S6.21 RESPONSE. It neither repairs the
delta=0 action nor changes the physical source metric. The proposed gate
is `P8-S6.23.PREPARED`; its report remains unpromoted until independent
review and read-only replay.

The actual parameters are `M,tau>0`, `0<delta<=10^-9`,
`u=T/tau`, `c=2+delta`, and `K=(tau*kcom)^2` in `[0,4]`.
The physical metrics have positive lapses `N_g=1,N_f=c` and scale factors
`a=(1+u^2)^2,b=2/(1+u^2)^2`. The signature/curvature convention is
`+---`, `R_B=-6(DH+2H^2)`, Einstein action `-M^2 R_B/2`.
The source is a linear external conserved transverse-traceless stress on
the actual `g` metric. No nonlinear source realization is assumed.

The theorem concerns the two-dimensional **source-free prepared sector**
on `|u|<=1/100`. It is selected by joint analyticity of `l` and
`q=Q/D`, where `D=(2+delta)(1+u^2)^4-2`, and by the even/odd canonical
initial light data `(l,l')=(1,0),(0,1)` at `u=0`. Existence is proved in a
weighted coefficient-l1 space on `|u|<=1/20,|delta|<=1/400`, uniformly for
complex `|K|<=4`. Complex and negative delta are only an analytic proof
device for the pole-cleared equations. They are not actual physical
parents. Joint analyticity is an explicit preparation condition, not a
condition on arbitrary incoming states or physical sources.

The strongest conclusions are as follows.

1. The full coupled tensor equations admit this analytic two-mode sector.
   A convergent bounded inverse and a rational contraction factor below
   `1/2` prove existence, not just formal inner polynomials. The norms obey
   `||l_e-1||<10v, ||q_e||<3v, ||l_o-u||<R/40, ||q_o||<R/5`,
   where `R=1/20,v=R^2`.
2. On the actual real slab the conserved restricted canonical symplectic
   form has `|Omega-1|<=3delta^3/5`; the physical g Wronskian satisfies
   `3/5<W_g<2`. The exact two-dimensional closed g equation has positive
   kinetic normalization `9/20<K_eff=Omega/W_g<11/6`.
3. This equation is analytic in `K`, including `K=0`, at which its
   zero-derivative coefficient vanishes exactly. At the center,
   `|G(0)/K-1-(14/33)delta|<=48001delta^2`, with the ratio at `K=0`
   defined by analytic continuation. This is not the locked value
   `1+(4/5)delta+O(delta^2)` and is not a cone or causality claim.
4. In the actual endpoint normalization of S6.21, for `1<=K<=4`, the
   prepared inclusion differs from the delta=0 regular-light inclusion
   by at most `200delta`. Fixed delta=0 light Cauchy data therefore have
   full fixed-slice output error at most `8600delta` times their norm.
   This is an upper rate, not a sharp leakage coefficient. A fixed old
   g-source adds `126000000delta` times its explicit source L1 norm.
5. An explicitly given delta-dependent, conserved g-only source supported
   strictly in `(-1/50,-1/100)` prepares each analytic solution exactly.
   Its source norm and O(delta) retuning have finite rational bounds for
   a specified flat cutoff. Its duration is fixed; it is not band-limited.

The physical source-free equation is a representation of a selected
solution subspace at each spatial momentum. It is **not** an off-shell
local EFT, a proof of arbitrary-source heavy-mode elimination, scalar or
vector health, a low temporal band, positivity, or an original C/D match.
In particular, the center coefficient requires the explicitly tuned
analytic preparation. A fixed source with O(delta) endpoint errors need
not enjoy the same center expansion. Neither P8(b) nor the full original
P8 problem is closed by this calculation.
