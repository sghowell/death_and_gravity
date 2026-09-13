# Full scalar metric, Hamiltonian and canonical gauge bridge

Write \(v=a^3,\ \tau_D=\operatorname{tr}Q_D\). Direct differentiation of the entire ADM metric gives
\[
 h_D=\begin{pmatrix}2n_D&-a^2\beta_D^t\\-a^2\beta_D&-a^2Q_D\end{pmatrix}.
\]
Its complete mixed chart derivative is
\[
 q_{00}=2n_Dn_G-2a^2\beta_D\cdot\beta_G,\quad
 q_{0i}=-a^2(Q_D\beta_G+Q_G\beta_D)_i,\quad
 q_{ij}=-\frac{a^2}{2}(Q_DQ_G+Q_GQ_D)_{ij}.
\]
The noncommuting order and all lapse/shift terms remain.

## Full scalar Hamiltonian vertices

For the five fixed reference features
\(F=(\Pi/v,\nabla\phi/a,m\phi)\), the first Hamiltonian vertex is
\[
 M_D=\begin{pmatrix}
 n_D-\tau_D/2&a\beta_D^t&0\\
 a\beta_D&(n_D+\tau_D/2)I-Q_D&0\\
 0&0&n_D+\tau_D/2
 \end{pmatrix}.
\]
Set \(b=(n_D\tau_G+n_G\tau_D)/2+\tau_D\tau_G/4\).
The full second vertex has momentum entry
\(\tau_D\tau_G/4-(n_D\tau_G+n_G\tau_D)/2\), mass entry \(b\), and gradient block
\[
 bI-(n_D+\tau_D/2)Q_G-(n_G+\tau_G/2)Q_D
       +(Q_DQ_G+Q_GQ_D)/2.
\]
All second shift entries vanish in these canonical ADM variables. This does not erase the nonzero shift terms in the metric chart \(q\); the two are different chain-rule objects.

Literal first and mixed derivatives of \(Ne^{-\operatorname{tr}Q/2}\), \(Ne^{\operatorname{tr}Q/2}e^{-Q}\), \(Ne^{\operatorname{tr}Q/2}\) and the shift term verify every entry. The quadratic matrix-exponential jet is exact for this derivative calculation, not a truncated finite-metric model. Independent finite differences use the entire numerical matrix exponential on noncommuting directions.

## Prepared gauge decomposition

For \(\xi=(\eta,\chi)\), direct evaluation of \(\mathcal L_\xi g\) gives
\[
 n_\xi=\eta',\quad \beta_\xi=\chi'-a^{-2}\nabla\eta,\quad
 Q_\xi=2H\eta I+\operatorname{symgrad}\chi.
\]
Let \(I_-f=\int_{-1/2}^t f\), and \(I_+f=-\int_t^{1/2}f\). Both derivatives equal \(f\).
For the source take
\[
 \eta_G=I_-n_G,\quad
 \chi_G=I_-(\beta_G+a^{-2}\nabla\eta_G),\quad
 G_{\rm syn}=(0,0,Q_G-2H\eta_GI-\operatorname{symgrad}\chi_G).
\]
Use \(I_+\) for the detector. No \(H\) or momentum is divided by. The source flow is identity near preparation and therefore carries the same initial scalar covariance.

## Complete scalar mode identity at nonzero transfer

Let \(\ell=k+p\). The background canonical evolution for \(z_k=(\phi_k,\Pi_k)\) is
\[
 z_k'=M_kz_k,\quad
 M_k=\begin{pmatrix}0&v^{-1}\\-v(m^2+|k|^2/a^2)&0\end{pmatrix}.
\]
For an arbitrary Fourier metric source, its complete perturbation from \(k\) to \(\ell\) is
\[
 \delta M_{\ell k}=
 \begin{pmatrix}
 i k\cdot\beta&(n-\tau/2)/v\\
 -v[(m^2+k\cdot\ell/a^2)(n+\tau/2)-\ell^tQk/a^2]&i\ell\cdot\beta
 \end{pmatrix}.
\]
The two shift momenta differ and both are necessary. The two-momentum Hermitian/symplectic identity is checked exactly, so the complete linear source preserves CCR.

A scalar gauge pullback gives
\[
 \delta\phi=\eta\phi'+\chi\cdot\nabla\phi,\qquad
 \delta\Pi=\eta\Pi'+\operatorname{div}(\chi\Pi)+a\nabla\eta\cdot\nabla\phi.
\]
Thus for the Fourier pair,
\[
 \delta z_\ell=B_{\ell k}z_k,\quad
 B_{\ell k}=\begin{pmatrix}
 ik\cdot\chi&\eta/v\\
 -v(m^2+k\cdot\ell/a^2)\eta&i\ell\cdot\chi
 \end{pmatrix}.
\]
Substituting the full gauge ADM source gives the exact identity
\[
 \delta M_{\ell k}=B_{\ell k}'+B_{\ell k}M_k-M_\ell B_{\ell k}.
\]
It includes the canonical-density gradient term. Since \(B=0\) on the common initial germ, uniqueness of the ordinary linear mode system identifies this with the retarded same-Cauchy-data tangent. Its full Gram tangent follows by the product rule, with both reverse blocks retained. This is a direct scalar verification, not an imported vector response formula.
