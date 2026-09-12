# P8 continuation: complete homogeneous trace anchor

S6.218 completes the homogeneous trace/trace Gaussian anchor, including the varying volume, longitudinal constraint, original dimensional finite matching and full metric contact. Original P8 remains OPEN. The [corrected tracefree current](assessment-2026-09-12-p8-affine-corrected-spatial-current.md) is unchanged; the [phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) still withdraws the old S207/S212/S213 physical endpoint identifications.

## Literal trace family and all-momentum comparison

For Q=(2/3)phi I and phi=epsilon Gamma, the effective scale is A=a exp(phi/3) and the volume is A^3. The constrained Hamiltonian has
\[
 K_T=A^{-1},\quad K_L=A^{-1}+k^2/(A^3m^2),\qquad
 V_T=Am^2+k^2/A,\quad V_L=Am^2.
\]
Both products are the same full frequency squared. The longitudinal temporal-constraint contribution is positive and is retained in both metric derivatives. This is a genuinely isotropic polarization decomposition, not a relabeling of the tracefree shear Hamiltonian.

With z=k^2/(A^2 omega^2), the normalized current vertices are diag((1-2z)/3,-1/3) and diag(1/3,-(1+2z)/3). Their first three amplitude bounds (1,1/9,2/27) lie below the preceding full-comparison displays. Exact mixed frequency and squeeze bounds check all inherited hypotheses, including the finite complement, nonzero initial reference mismatch, both covariance directions, full infinite tail and both marker branches.

The family has |epsilon|<=1/100, twelve bounded Gamma time jets and the unchanged common zero initial neighborhood. The actual-minus-fourth-order current has C0,C1,C2 bounds (1e77,1e94,1e111). These are full coordinate-density estimates; the volume is not held fixed.

Notation clarification: in the report's actual_M symbolic packet, A is the positive unperturbed base scale and ae=A exp(phi/3) is the effective scale. In the narrative Hamiltonian formulas, A denotes the effective scale. This is a naming distinction, not a change to frozen formulas or evidence.

## Original finite matching and homogeneous identification

Independent arbitrary-dimensional scalar Hamiltonian differentiation and rational radial moments retain the full A^d Jacobian. Weighted total derivatives are removed at general dimension, before the finite limit. The original prescribed physical local density is
\[
 \frac{A^3}{64\pi^2}
 \left[\frac52m^4+\frac53m^2R-\frac{R^2}{30}
 -\frac{Ric^2}{15}\right],
\]
modulo the compact four-dimensional Euler variation after the limit. No finite coefficient or preparation is changed.

Its complete scalar Euler current is A^3 L/(64 pi^2), where
\[
 L=\frac52m^4+10m^2H^2+\frac{20}3m^2H'
   -12H^2H'-8HH''-6(H')^2-\frac43H'''.
\]
The cosmological term gives a nonzero trace current and Hessian. All five local formal density Green identities hold; these are not symmetry of the retarded nonlocal kernel. All three local amplitude bounds are below1e11.

The complete determinant-component C0,C1,C2 displays are (1e78,1e95,2e111), with Taylor remainder at most epsilon^2 1e111. At the reference, the same covariance tangent, full second contact and finite-mode ODE uniqueness identify its first derivative with the actual homogeneous trace response. S176 supplies the actual vanishing source and first source variation there; this does not make the affine parent source free at finite amplitude.

Rotational covariance makes both ordered trace/tracefree cross anchors zero at P=0. It does not equate the two nonzero-transfer retarded cross kernels. Normalization to I/sqrt3 gives both scalar-coordinate conversion factors, with product3/4. The orthogonal full homogeneous lift satisfies
\[
 |H_0[D,G]|\le 10^{95}\|D\|_{L^2}Z_{130}[G].
\]
The canonical displays are3e-705 for unit trace and4e-705 for the full homogeneous lift. Thirteen source time derivatives remain; these weak bounds are not a reduced inverse or stability theorem.

## Independent checks and immutable fresh validation

Independent tests use literal ADM metric determinants, inverses, curl and temporal-constraint terms at zero, intermediate and large momentum. Separate finite covariance-flow derivatives retain the full first-response contact. Raw scale-factor Euler variations and integrated two-direction Hessians include the mixed second scale variation. General-dimension radial matching is independently derived.

Exploration first found a text-scope expectation mismatch, then a preflight exact-serialization issue: Python zero division had made twelve bound margins floats. Both were corrected before freezing without changing a physical coefficient. Independent science passed164 tests in13.89 seconds; integrated389 in14.55 seconds; exact preflight389 in3.91 seconds; final private389 in14.30 seconds; repository science389 in14.55 seconds. Lint passed.

Fresh original-SymPy ordinary replay passed414 tests in2405.23 seconds; CLI replay passed. The complete P8 regression passed43603 tests in4272.12 seconds with final exit code0. Its689-file snapshot SHA is
\[
 \texttt{fd7e54fe5e8479abd029d6f91718db2012c80e5f9f54b7316e967f40cb557790}.
\]
Only full regression used the audited exact-GCD adapter. Its128 original tuple comparisons passed, with final counters35851 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The51620-character native report was transferred in five checked chunks. Its20 fields contain56 named identities,110 scalar entries,46 gates,nine controls and167 rejected inputs. All19 frozen source/proof/test hashes and report SHA were verified:
\[
 \texttt{0c12575dca9e1ba1361343e7832e93e377e6198d7923ec49fe7d7f51e54b1153}.
\]
The exact23-file publication manifest includes those19 inputs, report, this audit, CLAIMS and README. Pending S219, private reduced-scalar work and unrelated P4/P9 changes are excluded. Continuous all-momentum, dimensional, dominated-C2 and Hilbert arguments are written proofs, not FORMALIZED.

## Remaining work

The nonzero-transfer scalar state/time/endpoint/contact and all-transfer response assembly is a separate successor under validation. A genuinely constrained reduced inverse, full nonlinear sourced-parent remainder, quantum-corrected background, physical cutoff/heavy-sector control and original V/G/B are not established here. Scoped P8(a), A.20–A.23 and the nine primitive statuses remain unchanged. There is no user-intervention blocker; continuation is active.
