# Entire actual response and restored original-frequency projection

All objects below use the same actual reference state, scalar prescription and fixed profile. Let \(\mathcal H_0\) be S6.244's full actual homogeneous spatial response. It includes the full Hamiltonian contact, finite heat action and fixed-profile variation.

For each external transfer \(P\), let \(E(P)\) be the complete actual-minus-computational-W6 ordered memory, \(F(P)\) the full fifth endpoint plus sixth bulk, and \(Q(P)\) the complete subtracted first-five endpoint remainder defined in the endpoint note. Let \(\mathcal L_{\rm UV}(P)\) be S6.243's full finite scalar transfer difference.

## Exact partition before estimates

At finite regulated integrals, subtract the zero-transfer equation from the equation at \(P\). The spatial Hamiltonian and fixed-profile contacts are transfer independent and cancel in this difference only. Decompose the entire ordered W6 time integral by the six-step identity. Subtract every unaveraged UV coefficient, retaining the complete remainders and restoring the same-scheme local difference.

The full dimensional dominated limit and frozen scalar coefficient match then give
\[
 \mathcal R(P)=\mathcal H_0+
       [E(P)-E(0)]+[F(P)-F(0)]+[Q(P)-Q(0)]
       +\mathcal L_{\rm UV}(P).
\]
Thus every nonzero transfer is anchored to an actual homogeneous response, not to a replacement WKB state. The anchor's contacts and profile have not been omitted. The exact algebraic repartition is checked independently in the packet.

The zero-transfer operator acts on the finite-dimensional symmetric tensor space and the time graph. It is not a claim that a nonzero spatially constant source is in \(L^2(\mathbb R^3)\). Its real Hilbert-space bound complexifies with the same norm and lifts componentwise to Fourier space.

For a prepared source, \(\sup_t|\partial_t^jG|\le\|\partial_t^{j+1}G\|_{L^2(I)}\) for \(j\le12\). This transfers the actual homogeneous source bound to \(H^{13}\) time. The other pieces require no more source derivatives than the final graph. Cauchy-Schwarz in transfer and Plancherel retain detector \(L^2\), with the factor \((1+P^2)^3\) wholly on the source.

## Full positive normalized sum

The anchor obeys
\(\|\mathcal H_0\|/\kappa<10^{45}n^2/\kappa\).
Use the exact \(C_Q,C_F,C_{\rm state}\) from the notes and the unchanged S6.243 local bound \(C_{\rm UV}<10^{-600}\). Then
\[
 \|\mathcal R\|/\kappa
 \le \frac{10^{45}n^2}{\kappa}
   +\frac{2C_Qn^2}{\kappa}
   +\frac{2C_F}{10^{98}\kappa}
   +\frac{10^{106}}{\kappa n^3}
   +C_{\rm UV}
 <\frac{10^{46}n^2}{\kappa}<10^{-350}.
\]
All differences are bounded by positive sums; no unproved cancellation is needed for this inequality. The finite local part appears exactly once.

## Definition and scope of the restored projection

For \(K\ge2m\), use the original
\[
 \chi_{P,K}(k)=1_{\Omega_k\le K}1_{\Omega_{P-k}\le K}.
\]
Let \(E_K,F_K,Q_K\) denote this mask applied to each COMPLETE convergent, unaveraged integrand. In particular \(Q_K\) includes the raw endpoint minus ALL its subtraction terms under the same mask. Let \(\mathcal H_{0,K}\) be the complete subtracted homogeneous projection from S6.244, with its finite local/profile pieces restored unchanged. Define
\[
 \mathcal R_K(P)=\mathcal H_{0,K}
   +[E_K(P)-E_K(0)]+[F_K(P)-F_K(0)]
   +[Q_K(P)-Q_K(0)]+\mathcal L_{\rm UV}(P).
\]
This explicit anchored, UV-subtracted projection is the regulator in the claim. It preserves the full finite local difference and profile. It is not an assertion that the bare projected Gaussian current alone equals this expression or obeys a covariant finite-\(K\) Ward identity.

At finite \(K\), the mask depends on angle and therefore an angular-zero UV term can contribute. No angular averaging that removes such a term is performed before masking. The previous note's nonzero hemisphere example is an explicit control.

## Full same-graph tail

The homogeneous tail is below \(10^{51}/(\kappa K^2)\). Both transfer differences retain their two complete tails. Since \(K>1,\ m<10^{99},\ R_K>K\),
\[
 \|\mathcal R-\mathcal R_K\|/\kappa
 \le\frac1K\left[
  \frac{10^{51}}{\kappa}
  +\frac{2C_{Q,K}n^2\,10^{99}}{\kappa}
  +\frac{2C_{F,K}}{\kappa}
  +\frac{10^{107}}{\kappa n^2\,10^{98}}\right]
 <\frac{10^{-250}}K.
\]
The exact positive constants, including the actual homogeneous bound rather than only its displayed upper bound, are serialized and checked. Large-small pairs, both removed legs, the unexpanded low region and every near logarithm remain.

## What this does not complete

This establishes a full actual FIRST spatial response at the fixed FLRW reference and a specified convergent projection. It does not establish a second derivative on finite inhomogeneous histories. For full ADM/common-clock assembly one must still retain the finite-cutoff local/shape bookkeeping and nonlinear contacts proportional to the nonzero projected mean tail. Those are not justified by taking this limit silently.

No quantum constraint inverse, stability or nonlinear response is inferred from graph smallness. Interacting light/mixed state and loop control, quantum gravitational decoupling, physical UV/Regge estimates and original V/G/B/P8 remain open.
