# Complete local current and actual homogeneous trace anchor

Let \(F(H,H')\) denote the physical finite polynomial in the previous note before \(A^3/(64\pi^2)\). Under \(\delta\phi=\chi\),
\[
 \delta A^3=\chi A^3,\qquad
 \delta H=\chi'/3,\qquad \delta H'=\chi''/3.
\]
Put \(D_H=\partial_t+3H\). The complete Euler current is
\[
 J_{\phi,\mathrm{local}}
 =\frac{A^3}{64\pi^2}\,L,\qquad
 L=F-\frac13D_HF_H+\frac13D_H^2F_{H'}.
\]
Explicitly, with \(H_j=\partial_t^jH\),
\[
 L=\frac52m^4+10m^2H_0^2+\frac{20}3m^2H_1
 -12H_0^2H_1-8H_0H_2-6H_1^2-\frac43H_3.
\]
The first term comes from varying volume. Keeping volume fixed would delete a genuine current and Hessian contribution.

At the reference, the local linear operator has coefficients, before \(a^3/(64\pi^2)\),
\[
 c_0=L,\qquad c_{j+1}=\tfrac13\partial_{H_j}L,\quad j=0,\ldots,3.
\]
For each \(k=0,\ldots,4\), exact calculation checks
\[
 c_k=\sum_{j=k}^4(-1)^j{j\choose k}D_H^{\,j-k}c_j.
\]
These are formal density Green identities for the local operator, not symmetry of a forward retarded nonlocal kernel. An independent Euler variation in raw scale-factor jets \(A,A',\ldots,A^{(4)}\), and a separately integrated two-direction full-chart Hessian, agree.

For a quantitative bound, take the absolute-coefficient polynomial of \(L\) at \(m=1000\). Evaluate its four \(H_j\) variables at the explicit base majorants plus \(1/300\). The volume is at most \(2(5/4)^6\), and \(1/(64\pi^2)<1/576\). Each amplitude derivative is bounded by applying
\[
 1+\frac13\sum_{j=0}^3\partial_{H_j}
\]
to this nonnegative polynomial; the1 is the derivative of the varying volume. The three exact raw bounds are
\[
 \frac{5625132891495009211}{169869312},\quad
 \frac{5625161941796211511}{169869312},\quad
 \frac{5625195992153673811}{169869312},
\]
each below \(10^{11}\). This requires only four source time jets for the local current and no time derivatives of the detector.

Adding the full actual-minus-fourth-order comparison bounds gives strict coordinate-density displays
\[
 (10^{77}+10^{11},10^{94}+10^{11},10^{111}+10^{11})
 <(10^{78},10^{95},2\,10^{111}).
\]
The family is \(C^2\) into \(C^0(I)\); its determinant-component Taylor remainder is at most \(\epsilon^2\,10^{111}\), zero at \(\epsilon=0\). This is not a full sourced-parent finite-amplitude remainder.

At \(P=0\) the S214/S216 full spatial first and second vertices coincide with the derivatives of this same constrained trace Hamiltonian. The covariance tangent includes both \(JM_GC\) and \(CM_GJ^t\), with the same zero initial tangent from the source germ. Finite-mode ODE uniqueness, the complete current comparison and original local matching identify the first derivative with the actual prescribed trace/trace response. No erroneous creation/annihilation endpoint bridge is used for this anchor.

At zero momentum the state, regulator and fixed local prescription are rotationally invariant. The spatial causal map therefore has the form \(\alpha Q+\beta\operatorname{tr}(Q)I\), where the coefficients are time operators. It sends trace to trace and tracefree to tracefree. Both ordered cross anchors vanish separately.

For \(T_0=I/\sqrt3\), write \(\phi=(\sqrt3/2)G_0\) and use the same detector conversion. The trace Hessian bound is multiplied by \(3/4\). Combine it with S194's tracefree anchor using the orthogonal Frobenius decomposition. The initial germ gives \(\sup|\partial_t^jG|\le\|\partial_t^{j+1}G\|_{L^2_t}\) for \(j\le12\). Real/imaginary Fourier components and the orthogonal tensor outputs square-sum, so there is no extra tensor-basis factor. Plancherel gives
\[
 |H_0[D,G]|\le10^{95}\|D\|_2Z_{130}[G].
\]
This is the lift of an explicit homogeneous kernel, not point evaluation of an arbitrary weak Fourier multiplier. The canonical factors \(4/\kappa\) give \(3\,10^{-705}\) for unit trace and \(4\,10^{-705}\) for the full lift, with the derivative loss retained.
