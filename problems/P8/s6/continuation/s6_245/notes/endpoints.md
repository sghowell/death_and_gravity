# Complete first-five endpoints in all momentum regions

Let \(Q(P)\) contain the first five complete computational-W6 endpoints. For \(r<m\), keep them unexpanded. For \(r\ge m\), subtract every unaveraged inverse-radius coefficient through degree \(4-j\) of endpoint order \(j=0,\ldots,4\). Each source jet is retained.

There are 15 order/degree slots, 35 source-jet entries and four ordered tensor geometric products. The exact coefficients are the 140 frozen S6.243 coefficients, as identified from all six full frequency iterates. A coefficient whose full angular average is zero is still part of the unaveraged subtraction.

Write \(U=m+P,\ \epsilon=1/100,\ L=2U/\epsilon=200U\). The full normalized row bound on \(|x|\le\epsilon/U\), after summing its source-jet coefficients, is
\[
 B_j=1000^2\left(\sum_rc_{jr}\right)R^{-j}.
\]
The five values are
\[
 (10^6,\ 6\cdot10^{10},\ 5.6\cdot10^{15},\
       7.28\cdot10^{20},\ 1.224\cdot10^{26}).
\]
A row means an operator on all source jets at that time; its scalar majorant is multiplied by their Euclidean norm. No source analyticity is assumed.

## Entire real row, not an unphased memory bound

The complete all-real endpoint of order \(j\) is bounded by
\[
 C^2\,2^{j+1}\Big(\sum_rc_{jr}\Big)R^{-j}
             \frac{\nu\mu}{(\nu+\mu)^{j+1}}.
\]
Since \(\nu+\mu\ge2m\), the first-five sum is bounded by
\[
 2C^2\sum_{j=0}^4\Big(\sum_rc_{jr}\Big)(mR)^{-j}\nu
 < A\nu,\qquad A=4\cdot10^{18}.
\]
The enormous physical mass is used only in this explicit positive inequality. This retains the phase-integration factors and all source jets.

## Far region \(r\ge L\)

The complete Cauchy tail beyond degree \(4-j\), evaluated at half the analytic radius or less, is bounded before radial measure by
\[
 2B_j\epsilon^{j-5}U^{5-j}r^{-4}.
\]
With measure \(r^2dr/(2\pi^2)\), integrating from the actual threshold \(L\) gives
\[
 \frac{B_jU^{4-j}}{18\epsilon^{4-j}}.
\]
There is no replacement of the transfer-dependent threshold by a fixed mass radius.

## Near region \(m\le r<L\)

Keep the full real row and every subtraction term separately. The entire raw contribution is bounded by
\[
 \frac{4AU^4}{9\epsilon^4}.
\]
For a Taylor coefficient of degree \(d\), if \(j+d<4\), its full radial integral is bounded by
\[
 \frac{B_j\,2^{4-j-d}U^{4-j}}
           {18(4-j-d)\epsilon^{4-j}}.
\]
If \(j+d=4\), retain the logarithm first and then use
\(\log(L/m)\le L/m\), obtaining
\[
 \frac{2B_j U^{5-j}}{18m\epsilon^{5-j}}.
\]
The explicit \(1/m\) in all five logarithmic slots is essential. It is not replaced by 1 before mass powers are counted.

## Unexpanded low region \(r<m\)

Here \(\nu<2m\), so the entire first-five endpoint integral is at most
\[
 A m^4/27.
\]
There is no normalized inverse-bare-momentum frame or small-momentum omission in this ball.

## Complete finite sum and source weight

For every \(0\le q\le6\), \(m\ge1\),
\[
 (m+P)^q\le8m^q(1+P^2)^3.
\]
For the logarithmic slots the retained inverse mass lowers the degree by one. Hence every full finite contribution has mass power at most four. Summing all five far terms, the near raw term, all 15 near Taylor terms and the unexpanded low term gives
\[
 |Q(P)|\le C_Q m^4(1+P^2)^3,\qquad
 C_Q=\frac{111268759571258000000000000000}{9}<10^{30}.
\]
Time/source-jet norms and detector norms are understood as above. Plancherel places the entire spatial weight on the source.

## Original two-leg removed union

Keep the actual comoving cutoff radius \(R_K=4\sqrt{K^2-m^2}\), not a new cutoff.

In the far region, \(r\ge200(m+P)\). If a pair is removed, then \(r>R_K/2\), including when only the second leg is removed. Integrating the same complete far envelope gives
\[
 \frac{2B_jU^{5-j}}{9\epsilon^{5-j}R_K}.
\]
In the near region, removal implies \(R_K<r+P<201U\); multiply each complete near bound, including every logarithm, by \(201U/R_K\). In the low region, removal implies \(R_K<r+P<U\); multiply the full low bound by \(U/R_K\).

The resulting maximum transfer power is six and maximum mass power five. The exact positive sum is
\[
 |Q-Q_K|\le C_{Q,K}\frac{m^5}{R_K}(1+P^2)^3,\qquad
 C_{Q,K}=\frac{67387527200968092800000000000000}{27}<10^{33}.
\]
Since \(R_K>K\) and \(m<10^{99}\), this gives the rational numerator used in the complete assembly.

## Why unaveraged UV terms must remain

Take mass-scaled \(r=6,P=8,R_K=10\), an admitted \(K/m>2\). The first leg is inside the ball, while the second-leg condition selects \(u\ge0\). Although
\[
 \int_{-1}^{1}\frac u2\,du=0,
 \qquad \int_0^1\frac u2\,du=\frac14.
\]
Consequently dropping an angular-zero row before applying the mask changes the projected integrand.

Independent quadrature retains the exact removed angular fraction and its breakpoints in 252 near/far/low cases, with \(P/m\) from 0 to 300 and \(K/m\) from 2 to 100000. All bounds hold; the largest near actual/bound ratio is approximately 0.387906. This checks region and cutoff bookkeeping but does not replace the full uniform proof.
