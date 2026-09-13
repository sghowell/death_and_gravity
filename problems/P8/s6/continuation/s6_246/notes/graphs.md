# Complete physical ADM and common-clock norm bounds

All constants below are new bounds for the actual heavy scalar. The old vector magnitudes and withdrawn response phases are not used.

## Synchronous maps on the one-sided graph

Both time primitives have \(L^2(I)\) norm at most 1. For \(j\ge1\), \(\partial_t^jIf=\partial_t^{j-1}f\). On complex radius \(1/4\) about the slab, \(|z|\le3/4,\ |1+z^2|\ge7/16\); hence
\[
 |H(z)|<7,\quad |a(z)^{-2}|<28,\quad
 |\delta(z)|<\tfrac12(16/7)^3<6.
\]
Cauchy gives complete derivative majorants \(b\,j!4^j\), with \(b=7,28,6\), respectively. Set
\[
 B_r(b)=\sum_{j=0}^{r}\binom rj b\,j!4^j.
\]

Every source derivative through thirteen is retained. The synchronous source's \(Z_{136}\) norm is bounded by
\[
 C_S=4[3+4B_{13}(7)+2B_{12}(28)]
 =62408287126207901212<10^{20}
\]
times \(U_{138}\). The factor 4 exceeds \(\sqrt{14}\); it bounds the complete Hilbert sum of time jets. The map adds at most two spatial derivatives to the lapse and one to shift, so \(H^8\) source controls the \(H^6\) synchronous spatial graph.

The detector synchronous \(L^2\) bound is
\[
 \|D_{\rm syn}\|_{L^2}\le[1+4\cdot7+2(1+28)]V_{02}(D)
 =87V_{02}(D)<100V_{02}(D).
\]
Because S6.245 needs no detector spatial derivative, two detector spatial derivatives suffice. In Fourier variables \(q=|P|^2\), use
\(q^2\le(1+q)^2\) for the detector and
\(q^2(1+q)^6\le(1+q)^8\) for the source. There is no detector time derivative.

## Actual scalar one-current constants

S6.240 proves \(|\rho|,|P|,|\rho'|,|P'|<B=10^{400}\). On this smaller slab,
\(a\le25/16,\ |H|\le8/5\). In particular
\[
 \|E\|_F^2/B^2\le(a_{\max}^3/2)^2+3(a_{\max}/2)^2<16,
\]
and
\[
 \|E'\|_F^2/B^2
 \le[a_{\max}^3(1+3H_{\max})/2]^2
       +3[a_{\max}(1+H_{\max})/2]^2<400.
\]
Thus \(\|E\|_F<4B,\ \|E'\|_F<20B\). The exact scale bounds, not just a loose \(a<2\), justify these density constants.

In the respective input norms, \(\eta\le1,\chi\le2\), so \(\xi\le3\). Its full spacetime derivative has norm below 4: separately bound \(\eta_t,\nabla\eta,\chi_t,\nabla\chi\) by \(1,1,2,2\). The gauge ADM tuple is below 13. These estimates require only two detector spatial derivatives and no detector time derivative.

The first chart has operator bound 6. The full second chart has bilinear bound 30, with component majorant \(2+8+6+6+4=26<30\). The complete source density-Lie term is bounded by
\((20+32+32)B=84B\). Including its chart term gives
\[
 C_{W_s}=6\cdot84+30\cdot4\cdot13=2064
\]
times \(B\).

For the synchronous source, \(Q_{\rm syn}\) and its spatial gradient are bounded by 13, and its time derivative by
\(1+16+8+4=29\), times \(U_{138}\). Therefore \(h_{G_{\rm syn}}\) and its spatial gradient are below 52, its time derivative below \(2\cdot4\cdot2\cdot13+4\cdot29=324\), and its complete spacetime derivative below 376. The entire detector Ward coefficient is
\[
 C_{W_d}=4(3\cdot376+2\cdot4\cdot52)+30\cdot4\cdot13^2=26456.
\]
Both ordered terms give \(28520B\), without deleting a chart or density contribution.

## Complete profile and ADM sum

The complete spatial profile is below \(8B\|D\|_F\|G\|_F\). The entire ADM profile is below \(64B\|D\|\|G\|\); for example \(a^3<8\) and \(|\operatorname{tr}Q|\le\sqrt3\|Q\|_F\) give the smaller coefficient \(8(2+2+3/4)<64\).

Let \(C_{245}\) be the exact normalized full spatial bound in the pinned report. The full physical response therefore has positive majorant
\[
 C_{\rm ADM}
 =100\cdot10^{20}(C_{245}+8B/\kappa)
       +(28520+64)B/\kappa<10^{-327}.
\]
The subtracted profile's magnitude is included before reconstruction, and the full restored profile is included afterward. No cancellation is needed for this inequality.

## Complete nonlinear-clock graph

The first map \(T\) has detector bound
\(1+\sqrt3<4\), since \(0<\delta\le1/2\).
The entire source bound is
\[
 C_T=4[3+4B_{13}(6)]
 =51511602072425569260<10^{20}.
\]
It does not increase either spatial degree or time-derivative degree. After the explicit full-reference cancellation of the two nonlinear mean contacts,
\[
 C_{\rm clock}=4\cdot10^{20}C_{\rm ADM}<10^{-300}.
\]
An independent Fraction calculation reconstructs the full positive sums and both exact Leibniz constants.

## Output-density convention

If a later force equation requires output normalization by \(a^{-3}\), it acts on the OUTPUT:
\[
 \langle D,(\kappa a^3)^{-1}R_{\rm clock}G\rangle
 =\kappa^{-1}R_{\rm clock}(a^{-3}D,G).
\]
Since \(a\ge1\), time multiplication by \(a^{-3}\) cannot increase \(V_{02}\). No commutation through the retarded source-time kernel is asserted. This observation preserves the weak bound, but does not give graph invariance, a feedback self-map or a quantum constraint inverse.
