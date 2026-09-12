# P8 continuation: corrected retarded phase and full ordered scalar UV input

S6.216 corrects the retarded annihilation-pair orientation and constructs the full arbitrary-trace spatial UV difference in the original covariant finite prescription. Original P8 remains OPEN. This is a local UV result, not completion of the scalar current or restoration of the old S213 assembled-current formula.

The [independent phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) remains in force. S207's physical branch identification, S212's old finite value and S213's old physical assembly are not re-endorsed by replaying their frozen artifacts.

## Correct physical branch

For annihilation amplitudes \(b_X=u^tM_Xu\), the sharp-detector/source-annihilation current has positive detector phase, negative source phase and positive imaginary part. Exact finite Fock matrices and the real covariance tangent fix this orientation independently of any UV counterterm. The creation-amplitude S198 integration identity remains correct in its own convention; the later annihilation-branch application was the error.

The corrected endpoint coefficient is \(i(-i)^j\), so each odd endpoint reverses relative to the frozen annihilation extraction. The sixth bulk requires a consistent phase correction as well. This checkpoint corrects the local UV coefficients; it does not infer the full bulk repair from a local finite adjustment.

## Full trace and longitudinal constraint

For \(Q\) with arbitrary trace \(\tau\), use \(B_Q=Q-\tau I/2\). The full temporal-constraint pair supplies a new longitudinal amplitude and its two ordered cross products and square. Ten field-strength geometries produce150 endpoint/source-jet products and350 retained radial coefficients.

All six spatial invariants are reconstructed at general dimension before taking the fixed-physical-source dimension derivative:
\[
 T=\operatorname{tr}(DG),\quad V=e^tDGe,\quad
 W=(e^tDe)(e^tGe),\quad S=\operatorname{tr}D\,\operatorname{tr}G,
\]
\[
 U_D=\operatorname{tr}D\,(e^tGe),\qquad
 U_G=\operatorname{tr}G\,(e^tDe).
\]
Three ordered scalar channels add105 entries. The general-dimensional reconstruction is nondegenerate, with determinant \(d^3(d-1)^3/2\). The two ordered cross channels are retained separately.

The full curvature comparison includes the exponential metric volume, scalar curvature, Ricci and Riemann squares and the continued Euler combination. With the original prescribed pole weights, every corrected physical logarithmic coefficient agrees with the full covariant Hessian. No finite coefficient was retuned to force that agreement.

Before correction, the three scalar discrepancies were exactly
\[
 2a''p^2G_0,\quad 2a'p^2G_1,\quad
 -2p^2(a''G_0+a'G_1)
\]
in unnormalized \(I/I,I/S,S/I\) channels, \(S=\operatorname{diag}(-1,-1,2)\). The calculated odd endpoint reversal accounts for all three. The original canonical tangent and independent curvature target were not changed.

## Original finite term and local norm

The original fixed radial prescription retains the lower comoving band, mode-dimension derivative, fixed-source invariant derivative and full continued curvature-volume contribution. In the normalization of the report its finite input is
\[
 -m^2F_2G/2+(1-\log2-\ell/2)F_4G
 -\frac{\partial_df_4|_3}{2\pi^2}
 +\frac{\partial_dH_{\rm pole}|_3}{32\pi^2}.
\]
The corrected physical pole has \(H_{\rm pole}=16f_4\). Pole agreement alone does not establish the finite term.

In particular, the corrected-minus-S212 tracefree finite value is
\[
 \delta F_{\rm finite}
 =\frac{p^2(3t^2+1)(-31T+30V)}{420\pi^2}G_0.
\]
The physical odd spatial coefficient vanishes, but its dimension derivative does not. Source-jet1 and2 differences are zero. The physical spatial pole and lower-band power term are unchanged.

The complete local arbitrary-spatial finite operator has the bound
\[
 |F_{\rm finite}[D,G]|\le10^6\|D\|_2Z_{24}[G],\qquad
 Z_{24}^2=\sum_{r=0}^2\|(1-\Delta)^2\partial_t^rG\|_2^2.
\]
Its exact weighted coefficient sum is \(66786194445073/464486400<10^6\). There is strict slack for nonzero norm product, with zero at zero norm product. The ordered coefficients satisfy the correct formal proper-time Green relation, including interchange of \(U_D,U_G\); this is not symmetry of forward retarded kernels.

## Independent checks and fresh replay

Independent tests include finite Fock/covariance phase identities, all six retarded integration steps, literal full ADM polarization pairs,120 integer-dimensional scalar field-strength entries,48 scalar curvature Hessians and8 noncommuting curvature cases. A separate full four-order WKB/Cauchy extraction checks all350 retained coefficients at three real dimensional probes. Negative controls retain the wrong-phase scalar discrepancies, nonzero evanescent finite change and omitted-volume errors.

An initial independent corrected run passed228 tests and failed one structural comparison of algebraically equal expressions. Before freezing, that assertion was changed to compare the exact simplified residual; no coefficient was changed. Expanded integrated science passed467 tests in99.37 seconds, preflight467 in36.50 seconds, final private467 in99.50 seconds and repository science467 in98.92 seconds. Lint passed before freezing.

Fresh original-SymPy ordinary replay passed492 tests in2494.19 seconds and CLI replay passed. The complete captured P8 regression passed42869 tests in4406.40 seconds with final exit code0. It retained685 test files; snapshot path-list SHA:
\[
 \texttt{d4ccd057685a6ccc4142abf33d560a9e184a00ecbd6c2c97d12e3723c67cdb45}.
\]
The full-only exact-GCD adapter passed128 original tuple comparisons. Final counters were35828 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct science, ordinary and CLI retained original SymPy.

The52237-character native report was transferred in five checked chunks. Its20 fields contain75 named identities,706 scalar entries,38 proof gates,nine controls and152 rejected inputs. All21 frozen source/proof/test hashes and the report SHA were verified:
\[
 \texttt{70e8538425cfab3035cea3c4c9d3a5da63797182625f79811bbc9ad2cd9ed815}.
\]
The exact25-file publication manifest excludes the pending S217 successor and unrelated P4/P9 work. No frozen predecessor scientific, proof, test or report bytes changed. Continuous dimensional, covariant and Sobolev arguments are written proofs, not FORMALIZED.

## Remaining work

Full corrected finite-endpoint/bulk/current reassembly and its original-regulator proof are separate from this local result. S215's known tracefree-input application is not restored here. The homogeneous trace anchor, complete scalar state/time/contact/tail assembly, genuinely reduced inverse, nonlinear quantum background and remaining original V/G/B obligations are not claimed complete by S216.

A corrected full-current successor is undergoing its own fresh replay. Scoped P8(a), A.20–A.23 and all nine original primitive statuses remain unchanged. No user-intervention blocker has been identified.
