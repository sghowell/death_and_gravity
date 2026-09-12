# Complete ordered Ward cancellation and the extra clock contact

Use S215's weight-one contravariant current density
\[
 E^{\mu\nu}=\operatorname{diag}(A,B,B,B),\quad
 A=-a^3\rho/2,\quad B=-aP_{\rm stress}/2.
\]
The letters A,B here denote current-density coefficients, not the fixed scalar retuning coefficients in the previous note. The actual one-point function is nonzero.

For xi=(eta,chi), the ADM gauge direction is
\[
 n_\xi=\eta',\quad \beta_\xi=\chi'-a^{-2}\nabla\eta,\quad
 Q_\xi=2H\eta I+\operatorname{symgrad}\chi.
\]
The full ADM second chart derivative includes q00=2nDnG-2a^2 betaD dot betaG and qij=-a^2(DG+GD)ij/2. It must be included with the current density.

For an arbitrary detector(nD,betaD,QD), directly contract the full Lie derivative of E and the complete chart. The ordered source Ward integrand reduces to
\[
 2n_D(A'\eta_G+A\,\operatorname{div}\chi_G)
 +2(a^2B+A)\beta_D\cdot\nabla\eta_G
 -a^2\operatorname{tr}Q_D
 [(B'+2HB)\eta_G+B\eta'_G+B\,\operatorname{div}\chi_G].
\]
The betaD dot chiG' terms cancel exactly between the density and chart pieces. So do the spatial QD:symgrad chiG terms. Setting the chart to zero leaves a nonzero inverse-square Fourier pole.

For a synchronous source(0,0,QsG), the separately derived detector conservation identity is
\[
 a^2B[\eta_D(\operatorname{tr}Q_{sG})'
      +\chi_D\cdot\nabla\operatorname{tr}Q_{sG}].
\]
Its H and spatial Jacobian terms cancel against the chart. Source covariance and detector conservation are distinct ordered identities; they do not equate two retarded cross kernels.

## Fourier scalar reconstruction

First work with smooth Fourier sources supported away from P=0. Put Pi=Phat Phat^T and, on each respective endpoint domain,
\[
 \eta=I n,\quad c=I b-|P|^2 I(a^{-2}\eta),\quad
 \chi=-iP c/|P|^2,\quad \beta=-iP b/|P|^2.
\]
Here Iminus integrates from the left for the source; Iplus=-integral_t^right for the detector. Both derivatives equal the integrand. Then c'=b-|P|^2 a^-2 eta and the original lapse/shift gauge directions are recovered exactly. The synchronous spatial output is
\[
 Q_s=2(\zeta-H\eta)I-2c\Pi,\quad \zeta=v+\delta n.
\]
Equivalently, after the coefficient-sector b=pv/2 substitution,
\[
 Q_s=2(v+\delta n-H I n)I
 -\Pi I p_v+2PP^T I(a^{-2}I n).
\]
No inverse magnitude of P occurs in this output.

Pair detector momentum-P with source momentum+P. The scalar source Ward form is
\[
 W_s=2n_D(A'\eta_G+A c_G)-2(a^2B+A)b_D\eta_G
 -6a^2\zeta_D[(B'+2HB)\eta_G+B n_G+B c_G].
\]
Set T_G=6(zetaG-H etaG)-2cG. The detector form is
\[
 W_d=a^2B[\eta_D T'_G-c_D T_G].
\]
Both opposite Fourier signs are checked directly against the full tensor identities. If the source chart is deleted, the coefficient of |P|^-2 is
\[
 2a^2 A b_D c'_G,
\]
which is not zero. An independent coordinate-density push calculation observes this divergence at momentum1e-15 and its cancellation in the full expression. A separate covariant coordinate pull checks the detector identity.

The complete ordered response is Rspatial(QsD,QsG)+Ws+Wd. The source gauge field vanishes on its unchanged initial germ. The detector gauge field vanishes on its final germ; at the initial surface the full differentiated flux vanishes because the source and current tangent do. The synchronous detector may be cut off inside that zero-output neighborhood without a time derivative in M. No state is reset.

## Additional clock-to-metric chart

The physical logarithmic spatial chart is
\[
 Q_{\rm physical}=2vI-\tfrac12\log R(N)I,\qquad
 R(N)=1+2\delta(N^{-2}-1).
\]
Its first lapse jet is2delta I and its mixed second lapse jet is(8delta^2-6delta)nDnG I. Thus, beyond the already included full ADM chart, the scalar clock Hessian requires
\[
 W'_{\rm clock}=
 \frac{3a^3P_{\rm stress}}2(8\delta^2-6\delta)n_Dn_G.
\]
Independent full metric composition detects this term at both endpoints and the bounce. The fixed coefficient one-point cancellation cannot be used to delete this Gaussian chain-rule term.

These identities use the original covariantly renormalized current and same prepared state. A sharp computational band is not diffeomorphism invariant. No finite-band Ward identity is used or asserted.
