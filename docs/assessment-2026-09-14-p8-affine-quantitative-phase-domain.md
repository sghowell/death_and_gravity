# P8 S6.269 assessment: evaluated finite nonlinear phase domain and same-seed volume comparison

S6.269 evaluates an explicit bounce-slice domain for the finite coherent
quantum construction. It keeps the same original pure nonzero-mode
reference, full canonical normalization, entire nonlinear spatial
reconstruction and all twelve source-pinned auxiliary invariants.
Original P8 remains OPEN.

Fix the2pi torus L=1 and the six wavevectors +-10^64 e_i along all
three axes. All eight physical configuration channels are retained:
48 canonical pairs and96 phase coordinates. Homogeneous variables
are external reference parameters, not an additional quantum state.
This finite family uses the original mode prescriptions but is not
literally the R3 state or a physical ultraviolet cutoff.

Let V be the full original pure unit-CCR covariance at u=0 and
z=S w with V=SS^T/2. For every complete physical row l,
||lS||^2=2lVl^T. This gives deterministic field bounds on
||w||<=2*10^20 without replacing V by a diagonal covariance or
re-minimizing the state. Both scalar boundaries, the central
canonical swap, both tensor polarizations, actual H SLE and the
original all-order three-mode Proca graph remain.

The full S264 scalar/tensor bounds, S242 actual-SLE comparison
and S190 matrix graph give explicit rational physical field
ceilings, including all six real Fourier basis functions,
longitudinal vector factors and kappa=10^800. The actual heavy
frequency lies in[10^98,10^100]. No floating representation of
kappa or its reciprocal is used.

A notation clarification is required for the frozen reference note:
Omega=sqrt(n+P^2/16) is the fixed S240 energy-balancing scale
Omega_star, not the instantaneous KG dispersion. The physical
omega(u)^2=n+P^2/a(u)^2, and the chi equation additionally retains
-3H'/2-9H^2/4. At the bounce its physical-field dispersion is
n+P^2. The S269 field proof correctly uses the S240 variables
(sqrt(Omega_star)chi,chi'/sqrt(Omega_star)) and their full mode
bound; it does not substitute Omega_star into the KG equation.
Thus the evaluated field and invariant bounds are unchanged.
The word actual beside Omega in that note refers to its value
at the fixed parameters and must not be read as a new dispersion
relation. The frozen note/report remain archived unchanged.

The nonlinear shape proof is applied explicitly in physical
Fourier A2, not by pretending this high-momentum datum is small
in the older A8 ball. All nonnegative submultiplicative Fourier
weights support the same full determinant contraction.
Here ||tau||A2<10^-274, ||v||A2<10^-329 and
||Q-I||A2<=2*10^-274. The complete metric and its inverse
differ fromI in A2 by less than10^-272.

The full variable-coefficient adjoint has the term
-(partial_l Q_ij)(partial_j lambda_i)/3. Its generic
off-slice divergence contacts are retained before imposing
divQ=0. The full mean-zero Neumann inverse has A0-to-A1
norm below4; constants are not inverted. The entire original
metric/vector/matter source gives ||D0||A0<10^-194.
The complete cotangent lift has metric momentum below10^-192
and tracefree momentum below10^-190.

The generic exact identity tr(gamma DQ*[B])=0 preserves the
trace momentum Pi_v. It does not erase the shear correction.
Thus the trace density is p=Pi_v/(3V), with the original
density factor, while the full shear includes the entire lift.
The matter shift keeps Pi_M=1/10 at the canonical origin.

All twelve actual invariants, including shear, electric,
magnetic, vector mass, both matter gradients and full
intrinsic curvature, are uniformly below10^-260 on the
ENTIRE support ball. This is strictly inside S266's
10^-250 box. The full positive lapse and temporal branches
therefore apply, without assigning Gaussian laws to
nonlinear density invariants. The mean-zero spatial
constraints are solved; residual translations are retained
as quantum constraints, not set to zero on each displacement.

Using every full |N_i|<2 bound along the convex invariant
segment gives |Nstar(z)-Nstar(0)|<=24*10^-260.
The original center root is within10^-375 of1.
The whole-strip physical derivative bound is |U_N|<=12,
not just the bounce-center value3/2. The full spatially
averaged exp(3v)R_full(0,Nstar)^(-3/4) then differs from
its actual center by less than10^-256.

For any smooth invariant cutoff equal1 on the whitened
core of radius10^20 and supported strictly inside
radius2*10^20, the unchanged named convex extension
Fext=F0+chi(F-F0) obeys ||Fext-1||infinity<10^-255.
F0 is its actual source-pinned center, not reset to1.
The positive unital coherent map therefore gives the
actual infinite-Hilbert operator bound
||Q_V(Fext)-I||<10^-255, also on the common
zero-translation-charge sector.

In the SAME seed, coherent phase probability has
covariance2V, so its whitened law is the standard
Gaussian on96 real coordinates. For t=R^2/2=5*10^39
the exact outside-core probability is
exp(-t) sum_(k=0)^47 t^k/k!. It is strictly positive.
An exact gamma Chernoff estimate and rational
exponential-series bound give tail<exp(-10^39).
This is a POVM concentration estimate, not compact
state support or a joint sharp law of noncommuting q,p.

For two admissible cutoffs with the SAME core, full F
and actual F0, the seed vector difference is less
than2*10^-255 exp(-10^39/2) and the seed mean
difference less than2*10^-255 exp(-10^39).
These are evaluated comparisons of two defined
finite bounce-slice volume regulators.

This is the UNCALIBRATED positive coherent observable.
It does not evaluate the distinct first-Weyl-calibrated
volume or Hamiltonian error, nor assert its positivity.
Every cutoff derivative and implicit physical Hessian
contact remains required for that separate calculation.
No evolved-state leakage, nonzero time interval,
original interacting mean, uniform continuum limit,
physical matching or omitted-loop estimate is inferred.

The initial field-norm probe passed in0.0018009589985013008
seconds and the full geometry probe in0.24556304200086743
seconds. The complete packet probe passed: geometry
14 checks/86 scalar entries/17 gates in0.27571312501095235
seconds; reference11/67/24 in85.3295963749988 seconds;
source11/11/7 in217.96183820802253 seconds; invariants
8/8/23 in1.4293607500148937 seconds.
The combined zero exit was consumed at2026-09-14 18:11:21 UTC.

Independent diagnostics passed218 science tests with
four complete proof packets deselected in2.10 seconds.
They include three increasing Galerkin ghost compressions
at three amplitudes, the wrong-adjoint counterexample,
12 generic non-diagonal momentum/trace tests, eight
correlated covariance fixtures, eight radial cutoff
integrals and the actual160-digit arbitrary-exponent
gamma-tail check. The finite matrices are diagnostics,
not full diffeomorphism or CCR regulators.

The contract has44 named exact checks,172 scalar
entries,76 proof gates, nine controls and167 rejected
inputs, with nine primitive and125 matching rows.
Matrix entries are not independent theorems. The
written Banach/operator arguments are not FORMALIZED.
Ruff fixed one unused import, formatted seven files
and passed final check/format check on all nine
Python files. All seven sibling export contracts passed.

The final private original-SymPy preflight passed in331.7834723331034 seconds, followed by all222 science tests in0.61 seconds with proof packets cached. RuntimeWarning was an error. The zero exit was consumed at2026-09-14 18:39:06 UTC. All18 sources,20 AST-derived fields, exact residuals, proof gates, controls, parent hashes and frontier rows passed.

All18 scientific files,93520 characters, were frozen at2026-09-14 18:40:02 UTC after matching every successful private-preflight hash. Their repository copies preserve those exact bytes.

Fresh repository original-SymPy preflight passed in333.14143333304673 seconds, followed by all222 science tests in0.84 seconds with proof packets cached. RuntimeWarning was an error. Its zero exit was consumed at2026-09-14 18:45:43 UTC. Every successful private/preflight/repository source hash matched.

The original-SymPy native rebuild protected5495 scientific inputs before and after construction. Its170331-byte ASCII report was transferred in15 exact chunks with14 observed CONTINUE acknowledgments, ending at offset168000 with2331 characters. Raw un-reserialized SHA256 isf4bc0f6843458556163553c192d89003388362ab22e863a9cbfd7fd542065107. Raw validation passed at2026-09-14 18:49:57 UTC for all18 native/private/repository source hashes,20 AST-derived fields, exact counts and the complete retained frontier. No native chunk was truncated; an overbroad transport diagnostic initially matched the scientific word untruncated, and the exact transport-marker, length and raw-SHA checks confirmed the intact report.

The ordinary original-SymPy package replay completed with exit status zero: all247 tests passed in3639.60 seconds (1 hour39 seconds). Both captured test files were present and unchanged. Its zero exit was consumed at2026-09-14 19:52:40 UTC.

The separate original-SymPy CLI replay completed with exit status zero, consumed at2026-09-14 19:52:40 UTC.

The separately captured complete P8 regression completed with exit status zero: all63451 tests passed in5478.78 seconds (1 hour31 minutes18 seconds). Its791-file snapshot SHA256 is6c881babf7367a77c7128a73b54547522ccd6a7f671e8bf397f9152e69bcd8e6, with613 static namespace ancestors and128 original-versus-adapter selfchecks. All captured files were present and unchanged. Final exact-GCD counters were68058 domain fallbacks,7388 exact descents and94 mixed fallbacks. Its zero exit was consumed at2026-09-14 20:22:35 UTC. S270 and later work were not in this captured snapshot. The full retained output has no actual transport-truncation markers.

Native/direct/ordinary/CLI use original SymPy. Only
the separately captured full regression uses the
audited exact-GCD adapter and original-path selfchecks.
Exact raw reports are not reserialized. Publication
requires every observed zero exit, exact source and
report hashes, scoped staging and independent remote
verification. Unrelated P4/P9 changes are excluded.

The [explicit finite quantum construction](assessment-2026-09-14-p8-affine-local-quantum-regulator.md),
[full spatial reduction](assessment-2026-09-14-p8-affine-nonlinear-spatial-reduction.md),
[quantitative auxiliary branch](assessment-2026-09-14-p8-affine-nonlinear-lapse-branch.md)
and [reference integrability boundary](assessment-2026-09-14-p8-affine-nonlinear-reference-volume.md)
retain their distinct scopes. S261's physical-parent binding
remains explicitly refuted and its frozen archive unchanged.
Original V/G/B/P8 remain OPEN; completed scoped P8(a) is unchanged.
Missing physical estimates remain research work, not a newly
asserted user-intervention blocker.
