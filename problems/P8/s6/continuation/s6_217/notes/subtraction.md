# Original mask, odd logarithmic shell and corrected finite term

It is not enough to claim that a sign preserves a norm: a sign could change a nondecaying subtraction. Here the full physical UV parity calculation checks exactly which terms change.

At \(d=3\), every unaveraged odd-endpoint product of nonlogarithmic radial order vanishes. The only seven nonzero odd products all have endpoint \(j=1\), radial degree \(3\), and no transfer dependence. In the geometry/source-jet notation of S212 their legacy values are

| Geometry | Source jet | Coefficient |
| --- | ---: | --- |
| 00 | 0 | \(m^2(1+t^2)^4(11t^2+1)/4\) |
| 01 | 1 | \(-m^2t(1+t^2)^5/4\) |
| 10 | 0 | \(m^2(1+t^2)^4(11t^2+1)/4\) |
| 10 | 1 | \(+m^2t(1+t^2)^5/4\) |
| TL | 0 | \(m^2(1+t^2)^4(3t^2+1)/2\) |
| LT | 0 | \(m^2(1+t^2)^4(3t^2+1)/2\) |
| LL | 0 | \(m^2(1+t^2)^4(11t^2+1)/2\) |

The two source-jet1 magnetic entries must not be dropped before the complete geometry sum. Full angular contraction gives the same source-jet0 coefficient
\[
 m^2a^2(53t^2+7)/15
\]
in each normalized tensor, vector and scalar tracefree channel. Independently, all105 physical tracefree spatial-difference UV coefficients are unchanged.

For \(|P|\le K/4\), the part lost by intersecting the original two balls lies at \(r\ge K-|P|\). A logarithmic radial term obeys
\[
 \log\frac K{K-|P|}
 \le\frac{|P|}{K-|P|}\le\frac{4|P|}{3K}.
\]
The first inequality follows by differentiating \(x/(1-x)+\log(1-x)\), which is zero at the origin and has nonnegative derivative \(x/(1-x)^2\). Thus the changed odd logarithmic shell contributes a decaying error, not \(K^3,K^2,K,K^0\). For large transfer use S209's separate whole-band bound; do not extend the small-transfer shell expansion outside its domain.

S209's conversion error is bounded by an absolute sum over individual endpoint/radial rows, with the same low/high split and S208 row majorants. Multiplying each row and its Taylor remainder by a unit sign preserves that sum, including the changed odd logarithmic shell. The complete all-transfer bound remains below \(10^{40}M Y/K\). All nondecaying terms come from unchanged even rows, so the actual S210 \(A_3,A_1\), complete \(K^2,K^0\) cancellations and physical spatial \(F_2,F_4\) remain valid.

Dimension differentiation is different: S216 finds an evanescent odd spatial coefficient. The required corrected finite term differs from S212 by
\[
 \delta F_{\mathrm{finite}}
 =\frac{p^2(3t^2+1)(-31T+30V)}{420\pi^2}\,G_0 .
\]
It is nonzero and cannot be deleted because the physical pole was unchanged. It comes from the original dimensional prescription, not an added finite counterterm. The source-jet1 and2 differences are zero.

Consequently the corrected original-mask subtractor is
\[
 K^3A_3+KA_1+F_2(K^2-m^2)/2+F_4\log(K/m)-F_{\mathrm{finite}}^+ .
\]
The fixed lower comoving band is retained. The odd shell remains in the quantitative error; it is not absorbed into an arbitrary finite constant.
