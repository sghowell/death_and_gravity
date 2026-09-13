# Complete actual-state-minus-computational-W6 difference

The physical state remains S6.240's single actual SLE. Let \(S\) be the exact KG comparison initialized with the complete W6 Cauchy data, and \(T\) the actual SLE. The instantaneous computational W6 family need not solve KG. The following comparison is algebra on its full feature values, not a claim that it defines a physical metric response.

## Full error and products

The exact-comparison/W6 Bogoliubov error is below \(10^{62}\Omega^{-13}\); the complete normalized W6 Cauchy vector is below 20. Its full normalized-vector error is therefore below \(4\cdot10^{63}\Omega^{-13}\). The five-feature readout has norm below \(10\sqrt\Omega\), giving error below \(4\cdot10^{64}\Omega^{-25/2}\).

Combine this with S6.242's actual/comparison estimate, with \(B=10^{83}\):
\[
 \|F_T-F_{W6}\|
 \le300B\Omega^{-21/2}+4\cdot10^{64}\Omega^{-25/2}
 <600B\Omega^{-21/2}.
\]
Retain \(\|F_T\|\le400\sqrt\Omega,\ \|F_{W6}\|\le200\sqrt\Omega\). The pair difference has coefficient
\(3\cdot600\cdot400=720000\), multiplying
\(B\sqrt{\nu\mu}(\nu^{-11}+\mu^{-11})\).
The full ordered four-feature difference is bounded by
\[
 720000(3\cdot400^2+3\cdot200^2)B
       \nu\mu(\nu^{-11}+\mu^{-11}).
\]
Both full volume factors increase this by at most \(64^2\); the resulting coefficient is below \(10^{20}B\).

The telescoping identity uses actual factors to the left of each difference and full reference factors to its right. It therefore includes all mixed terms, error squares and the four-error product. The reflected internal ordering is retained. No initial error or state-selection contribution is set to zero.

## All transfers with source-only spatial weight

The common positive frequency norm is Lipschitz in internal momentum. With
\(L=1+P/(4m)\), both \(\mu\le L\nu\) and \(\nu\le L\mu\) hold. Thus the preceding integrand is bounded by
\[
 10^{20}B L(\nu^{-9}+\mu^{-9}).
\]
With the original three-dimensional Fourier measure,
\(\int\Omega_k^{-9}\,d^3k/(2\pi)^3<n^{-3}\).
The \(P\)-minus-zero difference consequently has full memory bound
\(4\cdot10^{20}B L/n^3\), uniformly on the retarded time triangle.

Since \(L\le2(1+P^2)^3\), put this complete transfer weight on the source only. Cauchy-Schwarz in time on the unit slab and then in \(P\), with Plancherel, give the normalized bound
\[
 C_{\rm state}=\frac{10^{106}}{\kappa n^3}<10^{-1280}
\]
on detector \(L^2\) times source \(Z_{136}\). This argument does not transfer a symmetric detector/source norm from S6.242; it explicitly proves the one-sided graph used here.

## Original two-leg removed union

Let \(K\ge2m\), and remove a pair when either original frequency exceeds \(K\).

If \(P\le2K\), Lipschitz comparison gives
\(|\nu-\mu|\le P/4\le K/2\). On the removed union both frequencies exceed \(K/2\). The complete ninth-power one-leg tail is below \(43K^{-6}\). Retaining both legs, both transfers and the \(L\) weight gives a coefficient below \(400\cdot10^{20}B\). Use \(K^{-6}\le m^{-5}/K\).

If \(P>2K\ge4m\), retain all large-small pairs and use the full, untruncated estimate. Here \(L<P/(2m)\le P/2\), so \(L/P^2<1/(4K)\). The remaining source weight is stronger than required. The zero-transfer tail is retained separately in the same positive sum.

Both regions are bounded by the convenient common numerator
\[
 C_{{\rm state},K}
 =\frac{10^{107}}{\kappa n^2\,10^{98}}<10^{-1180},
\]
times \(1/K\). This is deliberately looser than the direct small-transfer \(m^{-5}\) estimate.

The full spatial Hamiltonian and fixed-profile contacts are transfer independent on the reference background. They cancel only in the \(P\)-minus-zero difference. They remain in the actual S6.244 homogeneous anchor, including its cutoff-tail accounting; they have not been dropped from the physical response.
