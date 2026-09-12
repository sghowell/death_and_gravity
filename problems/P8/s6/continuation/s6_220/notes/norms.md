# Infrared-safe weak scalar norms

For s=(n,zeta,b), use the unitary spatial Fourier transform and lambda=1+|P|^2:
\[
 V_{03}[s]^2=\int_I\int_{\mathbb R^3}\lambda^3|\hat s|^2,\qquad
 U_{138}[s]^2=\sum_{j=0}^{13}\int_I\int_{\mathbb R^3}
                 \lambda^8|\partial_t^j\hat s|^2.
\]
These are Euclidean component norms. The source has the original zero initial neighborhood; the detector has the original final zero neighborhood. Both time primitives have L2 norm at most1 on this unit slab. For j>=1, (I f)^(j)=f^(j-1).

For the detector, eta=I n and c=I b-|P|^2 I(a^-2 eta), with a^-2<=1, give ||c||<=lambda||s||. Since2sqrt3<4 and |H|<=2,
\[
 \|Q_s\|\le(14+2|P|^2)\|s\|
          \le16\lambda\|s\|.
\]
The original spatial current detector M has one spatial derivative, hence
\[
 M[Q_{sD}]\le20V_{03}[D].
\]
Pi has Frobenius norm and operator norm1. It does not incur an inverse radial momentum.

For source jets, the unchanged Cauchy bounds are H_j<=7j!4^j and (a^-2)_j<=28j!4^j. Write
\[
 B_r(b)=\sum_{j=0}^r{r\choose j}b\,j!4^j.
\]
Leibniz, the derivative identity for I and sqrt14<4 give
\[
 Z_{136}[Q_{sG}]
 \le4[6+4B_{13}(7)+2B_{12}(28)]\,U_{138}[G]
 <10^{20}U_{138}[G].
\]
The exact raw constant is S215's source raw constant plus12. The new argument applies after eliminating the inverse-momentum shift from the complete form, not by applying S215's unreduced beta norm to beta=-iP b/|P|^2.

## Both Ward terms and clock contact

The original physical stress and first derivative bounds are below1e30. On the slab a<2, |H|<=2 and |H'|<=4. Thus
\[
 |A|\le4e30,\quad |A'|\le28e30,\quad
 |B|\le e30,\quad |B'|\le3e30,
\]
and |a^2B+A|<=8e30, |B'+2HB|<=7e30. The explicit source expression in notes/ward.md contributes, in units1e30,
\[
 56+16+168+24+2(8+24)=328.
\]
Here ||eta||<=||s|| and the safe bound ||c||<=2lambda||s|| suffice.

From the same exact projector formula,
\[
 \|T_G\|\le22\lambda\|s_G\|,\qquad
 \|T_G'\|\le46\lambda\|s_G\|_{\text{time jets }0,1}.
\]
The latter includes H', H n and c'=b-|P|^2 a^-2 eta. No source endpoint term is dropped. The detector term contributes at most4(46+2*22)=360 in units1e30. Finally |8delta^2-6delta|<=5 and3a^3/2<12 give a clock contact below60e30.

Every product has at most two factors lambda, so Fourier/time Cauchy-Schwarz is controlled by the displayed V03 U138 weights. There is no detector time derivative. The full local correction is below748e30<1e35.

S219's full spatial response has bound2e95 M Z136. Consequently
\[
 |R_{\rm scalar}(D,G)|
 <[20\cdot10^{20}\cdot2e95+1e35]V_{03}U_{138}
 <1e117 V_{03}U_{138}.
\]
This includes the additional clock chart.

## Extension and regulator boundary

For smooth sources and detectors with a Fourier hole around the origin, all original metric and gauge directions are in the existing domain and the exact ordered identity applies. The estimate is independent of the hole radius. Such fields are dense in the stated R3 scalar Sobolev spaces; bounded bilinear continuation therefore defines the unique weak scalar extension. Pi can be assigned any value at the measure-zero origin. It is not direction-independently continuous there, nor does that assignment identify a global homogeneous solution.

This argument does not claim an unweighted L2 norm for beta itself, or construct a full finite-amplitude rough-metric/state family. It extends the reference bilinear form from its original smooth domain.

Replacing only the synchronous spatial response by S219's same homogeneous-anchored approximation, while retaining exact covariant Ward and clock terms, gives error
\[
 20\cdot1e20\cdot3e54\,V_{03}U_{138}/K
 <1e76 V_{03}U_{138}/K,\quad K\ge2000.
\]
This is a reconstructed spatial approximation. It is not identified with the literal sharp-band full lapse/shift current: that identification would require a finite-band Ward identity, which is not available.
