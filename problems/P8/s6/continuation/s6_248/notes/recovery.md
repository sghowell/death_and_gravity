# Both original constraints and the full current force subgraph

Let (E_n,E_v,E_b,E_sigma) denote the COMPLETE normalized original residuals, including current classical coefficients, pure Proca-plus-heavy Gaussian response and the applied force. With W=D+3H the full weighted adjoint is

E_eta=-W E_n+H E_v+W(delta E_v)+q E_b+ell_matter E_sigma,
E_w=E_v, E_c=-W E_b, E_r=E_sigma.

The variational boundary is a³[(E_n-delta E_v)eta+E_b c].
If every adapted residual vanishes, E_v=E_sigma0, then W E_b0. The ORIGINAL zero germ forces E_b0 by weighted first-order uniqueness, and then forces E_n0. All original lapse, shift and matter equations have therefore been recovered.

The germ cannot be removed. E_b=a^-3,
E_n=a^-3 integral_(-1/2)^t q(s)ds,
E_v=E_sigma0
give a nonzero adapted null residual with the actual background. They violate preparation. Equality only after the initial time is not the complete distributional boundary graph.

For a normalized physical three-force e=(e_n,e_zeta,e_b), the original clock tuple receives (e_n+delta e_zeta,e_zeta,e_b,0). Its adapted adjoint is
J_adapt e=(-W e_n+H e_zeta+q e_b,e_zeta,-W e_b,0).
The full equation is T_adapt x=-J_adapt e, with the same sign as the original forced saddle. A physical raw density is normalized at its OUTPUT time before this adjoint or any row primitive. The distinct source/detector clock maps and pure Gaussian response are exactly those in ward.md.

Recover the classical auxiliary Legendre variables
p_v=-6v'+6Theta n-3ell_matter sigma+2b,
p_sigma=sigma'+w_bg n.
They are not claimed to be the full quantum theory's canonical momenta. The exact current interface is

Z'=K_QG2 Z+F_QG2 g,
s=C_QG2 Z+D_QG2 g,
g=e+Qbar s.

The generic S222 saddle algebra is evaluated with the WHOLE QG2 coefficient map displayed in inverse.py: current Jc,A,Tcorr, actual Theta,E,ell_matter,delta,H and w_bg=-ell_matter E. S241 proves the actual positive lapse pivot and finite classical propagator hypotheses. No old QG1 numerical coefficient is silently reused.

The direct auxiliary matrix is
D=-u u^T/(2Jc)-(3/2)e_b e_b^T, u=(1,delta,0).
It has rank2 and is kept. With the prepared current classical propagator G0,
Z=G0 Fg and Taux=D+C G0 F, so

(I-Qbar Taux)g=e.

For every smooth prepared finite-ball e the coupled Euler inverse supplies such a g. Conversely any solution of this Schur equation reconstructs the full Euler solution, whose uniqueness has already been proved. This establishes both inverse identities on the stated smooth prepared finite-ball subgraph.

For a general mathematical phase drive h, let Z_h=G0h and use
(I-Qbar Taux)g=e+Qbar C Z_h,
Z=Z_h+G0Fg, s=C Z_h+Tauxg.
The full current Gaussian response and classical propagator preserve the smooth finite-ball class. This statement does not claim that every mathematical pair(e,h) is independently realizable by a physical preparation.

All resulting histories lie in the applicable derivative-losing weak response graph because every time derivative exists and bounded Fourier support gives every spatial Sobolev order. This proves an invariant SUBGRAPH and its bijection, not unrestricted completed-graph surjectivity, maximality, density in the graph norm or an all-Pmax norm bound. It does not infer a feedback self-map from an H13 weak estimate alone.

The coherent heavy mean is a separate exact block. In physical units its equation is
H_mean''+3Hubble H_mean'+(n+P²/a²)H_mean=j_H.
After dividing the action by kappa a³ its Euler operator is -kappa^-1 times this KG operator and its applied source is j_H/kappa. Thus no factor kappa is lost by setting the normalized heavy pivot to-1. S241's physical energy-scaled homogeneous propagator is below1000 for all P, including0; Duhamel must use the corresponding physical source scaling. Its decoupled coherent mean is consistent with, and does not eliminate, the nonzero heavy Gaussian metric response already included in the coupled four-block construction.
