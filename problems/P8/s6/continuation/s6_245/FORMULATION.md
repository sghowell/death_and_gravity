# Complete actual scalar spatial response, without changing the parent

## Fixed data and observable

Keep S6.240's actual heavy SLE, its single common preparation, physical CCR and complete covariant minimal-scalar prescription. The parameters are
\[
 n=10^{200}/512+2,\quad m=\sqrt n,\quad\kappa=10^{800},
 \quad a=(1+t^2)^2,\quad H=4t/(1+t^2).
\]
In particular \(10^{197}<n<10^{198}\) and \(10^{98}<m<10^{99}\). The physical clock slab is \(I=[-1/2,1/2]\). The fixed full reference profile is not re-minimized or varied.

The perturbation is the first spatial variation \(h=a^2\exp Q\), evaluated at \(Q=0,N=X=1\). Sources are smooth real symmetric tensors, initially spatially Schwartz, with a common zero germ at \(t=-1/2\). They need not vanish at the observation endpoint. Detectors are real \(L^2\) symmetric tensors. Closure extends the bilinear response to the prepared source graph stated in the README.

Use the original Fourier measure \(d^3P/(2\pi)^3\). Complex Fourier coefficients mean the linear complexification of the real tensor kernel. The physical sharp operation changes the internal mode branch, not external coefficients. An antilinear imaginary-part operation must not be applied to the external complex direction.

## Claim

The full actual linear spatial response, divided by \(\kappa\), is bounded by \(10^{-350}\) on detector \(L^2\) times source \(Z_{136}\), with all external transfers admitted. S6.244 supplies the actual homogeneous anchor, including its full contact and fixed profile. S6.243 supplies precisely the same-scheme complete scalar UV transfer difference. The computational full six-iterate WKB family is used only to decompose the actual kernel; it is not asserted to be a physical evolved state or a new metric-response prescription.

Let \(\Omega_k^2=n+|k|^2/16\) and
\[
 \chi_{P,K}(k)=1_{\Omega_k\le K}\,1_{\Omega_{P-k}\le K},\quad K\ge2m.
\]
The restored projection of [assembly](notes/assembly.md) uses this mask on every complete unaveraged convergent piece and retains all full finite local/profile contributions. Its same-graph tail is less than \(10^{-250}/K\).

The regulator is a mathematical projection used to estimate convergence. It is not a physical momentum cutoff. The statement does not identify the restored projection with a bare finite-cutoff Gaussian Ward system.

## Required full terms

The construction retains the complete actual-SLE-minus-W6 four-factor difference, all six integration-by-parts steps, the fifth endpoint and sixth bulk with their own signs, the unexpanded low band, the whole near-band raw endpoint, all 15 Taylor slots and 35 source-jet entries, and the far Cauchy remainder. Four ordered geometric products give 140 coefficients before angular simplification. No angular-zero coefficient may be removed before the two-leg mask.

The dimension limit continues the complete scalar frequency and geometric coefficients, holds all six physical invariants fixed and uses both mode branches at the same dimension. It does not assume a positive Hilbert norm in noninteger dimension.

## Explicit exclusions

No second derivative on a finite inhomogeneous history tube is established. Full prepared ADM/common-clock response, finite-cutoff Ward and local/shape bookkeeping, a compatible quantum constraint inverse, stability, nonlinear feedback, interacting light/mixed states and omitted loops, quantum gravitational decoupling, physical UV/Regge estimates and original V/G/B/P8 closure remain open.

A tiny graph norm is not an inverse theorem. The old finite inverse statements, old vector multiplicities and withdrawn odd-endpoint phases are not transferred. The S6.241 heat action is not added again.
