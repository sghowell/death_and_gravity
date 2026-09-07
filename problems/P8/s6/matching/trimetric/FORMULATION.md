# P8-S6.13.AUXILIARY — a chosen trimetric parent applicability screen

This new sibling pins S6.12 as context and leaves every old parent,
composite and CD/M1 certificate unchanged. It audits the specific
auxiliary model selected in arXiv:1804.04671, not every trimetric theory
or every completion of an effective matter coupling.

## Specified action and physical metric

Use real invertible vierbeine e,v,u with positive oriented determinants,
`g=e^T eta e`, `f=v^T eta v`, `h=u^T eta u`, `eta=diag(1,-1,-1,-1)`.
The action, in the repository curvature/sign convention, is

`S=-(G/2) integral sqrt|g|R[g] -(F/2) integral sqrt|f|R[f]`

`-2 integral det(u)[B+tr(u^-1 Q)] +epsilon S_m[h,psi]`,

`Q=p_g e+p_f v`, `G,F>0`, `B=beta0g+beta0f`.

There is no Einstein term for h, and no separate g/f cosmological term.
This is the source's choice `m_h²=0`, `beta_n^g=beta_n^f=0` for n>=2.
Links p_g,p_f are constant; the no-flat theorem requires at least one
nonzero link and does not require either sign or B!=0. Matter depends
only on the actual auxiliary vierbein/metric u,h. For the explicit
source dictionary use a positive canonical scalar
`L_m=det(u)[h^{mu nu} partial_mu psi partial_nu psi/2-V(psi)]`.

The parent square-root/symmetrization chart is additionally required for
the original metric interpretation. Ambient noncommuting matrix fixtures
audit variations; they are not asserted to solve those constraints.

| Quantity | Mass dimension |
| --- | --- |
| G,F | 2 |
| B,p_g,p_f,beta4g,beta4f,V | 4 |
| e,v,u,g,f,h,epsilon | 0 |
| psi, partial_mu psi | 1,2 |

Here epsilon is the paper's matter-strength parameter, **not** the
old asymmetric composite-family parameter 1/Y.

## Exact result and prerequisites

For constant nonsingular vierbeine the Einstein tensors vanish, while

`E_e=-2 p_g det(u) u^-T`, `E_v=-2 p_f det(u) u^-T`.

Thus no regular constant flat triple solves the full equations if either
link is nonzero. This does not divide by B or use the u equation. A
constant matter potential only replaces `B` by `B+epsilon V/2`; it
cannot cancel either dynamical-metric equation. The result rules out
the stipulated regular Lorentz-invariant Minkowski vacuum prerequisite
for this chosen parent. It is not a criticism of the source's ghost-free
degree-of-freedom construction, which need not possess that vacuum.

Where B and Q are invertible, the source-free auxiliary solution is
`u0=-3Q/B`. Its effective potential has restricted coefficients
`beta_n=C r^n`, `C=-54 p_g^4/B^3`, `r=p_f/p_g` when p_g!=0, and fixed
composite weight ratio b/a=r. It is not the old beta2-only parent.
The full physical matter metric is h, not its leading composite value
h0. First-order canonical-source and Lorentz-constraint corrections are
derived explicitly; no generic replacement by the frozen metric
composite action is certified.

## Deliberately separate countercontrol

Adding `-2 beta4g det(e)-2 beta4f det(v)` changes the parent. With
`p_g=p_f=q>0`, `B=-6q`, `beta4g=beta4f=-q`, zero scalar vacuum energy,
and `e=v=u=I`, every flat equation and the vacuum density vanish.
The relative flat quadratic Fierz–Pauli mass is

`m_FP²=q(1/G+1/F)>0`.

For equal Einstein coefficients this is `2q/G`, in the actual physical
clock h=eta. This separately named extension is a counterexample to
generalizing the no-vacuum result across added bare cosmological terms.
It is not an established rolling CD/M1 match, full-parent stability,
quantum completion, or cutoff hierarchy.

The [proof](notes/proof.md) distinguishes the exact no-vacuum statement,
formal source expansion, erroneous truncated cancellation, and extended
flat quadratic control. The [source audit](notes/sources.md) states the
literature and metric/vierbein boundaries. S6 and P8 remain open.
