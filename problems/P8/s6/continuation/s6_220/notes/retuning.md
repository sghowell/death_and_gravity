# Full fixed QG1 coefficient Hessian

The current candidate is CD-REG-AFFINE-ISO-QG1, not the old bare scalar coefficient. S182 defines the complete correction
\[
 \Delta F=-p_v+T(X)[-p+p_v-(r+p)(X-1)/2],
 \qquad T(X)=\frac{X^{1024}}{X^{1024}+(1-X)^{1024}}.
\]
The functions r=rho_ref/kappa and p=P_ref/kappa are evaluated once from the original complete reference history, state and prescription, and then fixed as scalar coefficient functions. They are not recomputed on each varied history. The constant p_v is the original fixed vacuum pressure normalization.

At the clock, the zeroth and first X jets are A=-p and B=-(r+p)/2; X jets2..1023 vanish. This is sufficient for a local second variation but is not a replacement of the full global coefficient. The literal clock lapse density is
\[
 N R(N)^{-3/4}[A+B(N^{-2}-1)],\qquad
 R(N)=1+2\delta(N^{-2}-1).
\]
Including the full spatial volume exp(3v), its linear term is
\[
 3Av+T_c n,\qquad T_c=(1+3\delta)A-2B=r-3\delta p,
\]
and its quadratic addition is exactly
\[
 \Delta J\,n^2+3T_cnv+\frac92Av^2,\qquad
 \Delta J=\frac{21\delta^2-3\delta}{2}A+(1-6\delta)B.
\]
Keeping only DeltaJ would delete genuine mixed and volume terms. Independent high-precision differentiation of the complete order1024 switch and metric volume checks all terms with nonzero calibration profiles. Those fixtures are algebra tests, not a change to the actual fixed profile.

On the unit slab, 32/125<=delta<=1/2 and the established fixed stress jets through5 are below epsilon=1e-770. The original lapse numerator has positive coefficients and constant1215, with denominator800(1+t^2)^18. Hence
\[
 |\Delta J|<4\epsilon,\quad |T_c|<(5/2)\epsilon,\quad
 J_c\ge\frac{1215}{800(5/4)^{18}}-4\epsilon>1/100.
\]
The profile is smooth, but no numerical bound for its time jets6..13 is inferred from that five-jet estimate.

The coefficient-sector one-point terms are not separately zero at the reference. They cancel the actual Gaussian one-point function in the retained combined equations. This cancellation is not permission to set the Gaussian one-point function to zero in a Hessian chart rule. The full Gaussian response and both chart contacts must still be added to the coupled linearization.
