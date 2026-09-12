# All six endpoints and the full bulk

For the positive detector phase define
\[
 K_+[f](t)=e^{i\Theta(t)}\int_{t_-}^t e^{-i\Theta(s)}f(s)\,ds,
 \qquad {\cal D}f=\partial_s(f/\Omega),\quad \Theta'=\Omega>0.
\]
With the unchanged zero initial source germ, integration by parts gives exactly
\[
 K_+[f]
 =\sum_{j=0}^5 i(-i)^j\,\Omega^{-1}{\cal D}^jf
   +(-i)^6K_+[{\cal D}^6f].
\]
For \(j=0\), use \(\partial_s e^{-i\Theta}=-i\Omega e^{-i\Theta}\); the recursive remainder is \(-iK_+[{\cal D}f]\). This proves the displayed formula by induction without discarding any upper endpoint. The phase convention gives \((-1)^j\) times each legacy annihilation endpoint after the appropriate imaginary part. In particular the fifth endpoint changes sign. The sixth bulk must use \(K_+\) and the correct positive-imaginary contraction; there is no universal rule multiplying the entire old oscillatory bulk by a real sign.

Partition the first five endpoints using the unchanged comoving split. If \(E_j\) is a legacy complete endpoint and \(U_j\) its high-band radial subtraction, define
\[
 E_j^+=(-1)^jE_j,\quad U_j^+=(-1)^jU_j,\quad
 Q^+=\sum_{j=0}^4\left[E_{j,\mathrm{low}}^+
                  +(E_j^+-U_j^+)_{\mathrm{high}}\right].
\]
The low band is not expanded. All finite cells and the original two-leg mask remain. Write \(F^+\) for the corrected fifth endpoint plus sixth bulk and \(R^+\) for the actual-state-minus-unit comparison in the same convention. Then
\[
 J_{\mathrm{actual},K}=Known_K^++C_{\mathrm{unit},K}+U_K^+,
 \qquad Known^+=R^++F^++Q^+.
\]
This is an exact repartition of the original finite-regulator current, not a fitted definition of a new current.

The S197/S201 comparison estimates are sums of absolute values of the two mixed pair differences and the error square. Consistent branch conjugation preserves them. The inverse-phase estimates for the fifth endpoint and sixth bulk are likewise absolute; at physical dimension both phases have unit modulus. Their combined bound and tail remain \(2\,10^{48}\) and \(2\,10^{52}/K\).

For \(Q^+\), each complete row and its subtraction have the same unit sign. Thus the S208 low, near, far and original removed-union estimates apply to their absolute values unchanged. Adding these established bounds gives the strict safe displays
\[
 |Known^+[D,G]|<5\,10^{48}M[D]Y[G],\qquad
 |Known_K^+-Known^+|<5\,10^{53}M[D]Y[G]/K .
\]
Use the same pointwise kernel estimates at \(P=0\) before Fourier Cauchy-Schwarz. This does not assume that point evaluation of an arbitrary weak multiplier is bounded.

Only rebase from \(R^++F^+\). Do not add the earlier Q203/Pfin pieces or the S204 local-action target again. The contact is retained in the homogeneous anchor. These are Gaussian weak-current estimates, not strong operator-norm stress differentiability.
