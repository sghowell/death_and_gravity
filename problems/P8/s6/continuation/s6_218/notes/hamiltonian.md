# Full trace Hamiltonian and current vertices

The independent S214 constrained ADM Hamiltonian has \(h=a^2e^Q\), \(\tau=\operatorname{tr}Q\), \(B_Q=Q-\tau I/2\). With \(Q=2\phi I/3\), \(\tau=2\phi\) and \(B_Q=-\phi I/3\). The electric and magnetic terms scale by \(e^{-\phi/3}\), the mass term by \(e^{+\phi/3}\), and the temporal-constraint term by \(e^{-\phi}\). None is a shear-volume cancellation.

For a fixed momentum direction the physical T/T/L Hamiltonian blocks are
\[
 K_T=A^{-1},\quad K_L=A^{-1}+\frac{k^2}{A^3m^2},\qquad
 V_T=Am^2+\frac{k^2}{A},\quad V_L=Am^2.
\]
Both products \(K_sV_s=\omega^2=m^2+k^2/A^2\) are exact. The positive longitudinal constraint is essential to the second equality. The transverse degeneracy is harmless: this isotropic family is genuinely diagonal in the fixed polarization projectors, not a diagonal approximation of a noncommuting shear family.

The balanced scalar frame is
\[
 T_s=\operatorname{diag}(\sqrt{K_s/\omega},\sqrt{\omega/K_s}).
\]
Its rotation vanishes. Put \(H_A=A'/A\), \(p=\omega'/\omega=-zH_A\), \(z=k^2/(A^2\omega^2)\in[0,1]\). Direct logarithmic differentiation gives
\[
 S_T=(H_A+p)/2,\qquad S_L=(H_A-p)/2.
\]
At zero momentum the two species have identical Hamiltonians, squeezes and current vertices; no singular choice of a longitudinal axis affects the answer.

For the scalar detector \(\chi\), the complete normalized current vertex \(G=T^tM_\phi T/\omega\), with \(\chi\) factored out, is
\[
 G_T=\operatorname{diag}((1-2z)/3,-1/3),\qquad
 G_L=\operatorname{diag}(1/3,-(1+2z)/3).
\]
These are the full physical vertices, not just leading high-frequency terms. Since \(z_\epsilon=-(2/3)\Gamma z(1-z)\), the only nonzero first derivatives have value \(4\Gamma z(1-z)/9\), and the second derivatives have value \(-8\Gamma^2z(1-z)(1-2z)/27\). Thus
\[
 \|G\|\le1,\quad \|G_\epsilon\|\le1/9,\quad
 \|G_{\epsilon\epsilon}\|\le2/27.
\]
All frame and frequency derivatives are included by differentiating these exact normalized expressions. They are bounded by the previous displays \((2,25,650)\).

The current's first source variation includes \(-\operatorname{tr}(M_{\phi\phi}\Gamma C+M_\phi C_\epsilon)/2\). The temporal part of \(K_{\phi\phi}\) is \(+k^2/(A^3m^2)\), not zero. Deleting it changes the normalized longitudinal vertex by the nonzero amount \(-z\). Full literal ADM matrices and finite covariance-flow derivatives independently check these statements.
