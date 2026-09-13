# Original projected mean and a precisely defined Ward-completed approximant

## Full actual mean, same original mask

Keep \(\Omega_k^2=m^2+|k|^2/16\), independent of time and the varied history. For \(K\ge2m\), define \(\rho_K,P_K\) by cutting the COMPLETE actual-SLE-minus-ad4 mode integrands at \(\Omega_k\le K\), then restoring the full unchanged finite local terms. This is not the bare divergent mode stress.

S6.240's complete integrand and its first time derivative are bounded by \(10^{301}\Omega^{-5}\). With \(r=4\sqrt{\Omega^2-m^2}\),
\[
 r^2\,dr=64\Omega\sqrt{\Omega^2-m^2}\,d\Omega
 \le64\Omega^2\,d\Omega.
\]
The original Fourier measure therefore gives, for both stress components and their first derivatives,
\[
 |\partial_t^j(\rho-\rho_K)|,\ 
 |\partial_t^j(P-P_K)|
 <\frac{16\cdot10^{301}}{9K^2}
 <\frac{10^{310}}{K^2},\qquad j=0,1.
\]
No mass-dependent or time-dependent replacement cutoff is introduced. The finite local terms cancel from this difference because they remain full on both sides.

## Conserved projected reference mean is not full cutoff covariance

For a physical exact scalar mode the full energy and pressure obey their pointwise mode Ward identity. For each complete adiabatic order \(2j=0,2,4\), write its general-spatial-dimension readouts as
\(a^{-d}\omega^{1-2j}e_j\) and \(a^{-d}\omega^{1-2j}p_j\), with \(\lambda=\omega'/\omega\). The exact packet rechecks all three identities
\[
 \partial_t e_j+(1-2j)\lambda e_j+dHp_j=0
\]
using the full scalar coefficients.

The fixed \(\Omega\) mask has no time derivative. Each subtracted mean integral and the full finite covariant stress consequently give
\[
 \rho_K'+3H(\rho_K+P_K)=0
\]
at the reference. This is a homogeneous MEAN identity. It does not make a sharp spatial projector covariant under arbitrary metric perturbations.

## Explicit ordered approximant

Let \(R_{245,K}\) be S6.245's anchored, complete unaveraged two-leg spatial projection, retaining all subtraction terms before masking. Let
\(E_K=\operatorname{diag}(-a^3\rho_K/2,-aP_KI/2)\).
Define
\[
 R_{{\rm ADM},K}(D,G)
 =[R_{245,K}-P_{QQ}](D_{\rm syn},G_{\rm syn})
       +W_s(D,\xi_G;E_K)+W_d(\xi_D,G_{\rm syn};E_K)
       +P_2(D,G).
\]
Both Ward terms use the same projected mean; the entire fixed profile remains full. The construction extends the specified synchronous block using its conserved reference mean and the ordered Ward formulas.

This is a Ward-completed computational approximant. It is NOT asserted to be the bare projected full scalar metric response. In particular no local/shape counteraction implementing that identification has been proved here. A prescribed convergent approximation must not be confused with a physical-cutoff bridge.

Under the unchanged nonlinear clock, retain the actual mean-tail contact:
\[
 R_{{\rm clock},K}(D,G)
 =R_{{\rm ADM},K}(TD,TG)
       +3a^3(P_K-P)(4\delta^2-3\delta)n_Dn_G.
\]
It is not permissible to substitute \(P_K=P\) just because the full reference mean matches the profile.

## Complete same-graph error

Let \(C_{245,K}<10^{-250}\) be the exact numerator of the S6.245 spatial tail, and let \(B_K=10^{310}\). Linearity of both ordered Ward terms and the full \(C^1\) mean bound give
\[
 |R_{\rm ADM}-R_{{\rm ADM},K}|/\kappa
 \le\left[
 100\cdot10^{20}\frac{C_{245,K}}K
       +\frac{28520B_K}{\kappa K^2}\right]V_{02}U_{138}.
\]
The identical fixed profile cancels from this comparison; neither mean correction is omitted.

The nonlinear contact has magnitude at most \(24B_K/(\kappa K^2)\). Indeed \(a^3<8\), and for \(0\le\delta\le1/2\),
\[
 0\le3\delta-4\delta^2\le9/16<1,\quad
 9/16-(3\delta-4\delta^2)=(2\delta-3/4)^2.
\]
After both complete clock maps, and using \(K^{-2}\le K^{-1}\),
\[
 C_{{\rm clock},K}
 =4\cdot10^{20}\left[
    100\cdot10^{20}C_{245,K}+\frac{28520B_K}{\kappa}\right]
       +\frac{24B_K}{\kappa}<10^{-200}.
\]
This proves the displayed \(1/K\) rate on the SAME \(V_{02}\times U_{138}\) graph.

All original spatial two-leg regions, including angular-zero UV terms before masking, remain in the pinned S6.245 block. The original one-leg mean tail and its full nonlinear clock contact are additional retained terms, not replacements for that spatial cutoff accounting.
