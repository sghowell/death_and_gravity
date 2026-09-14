# P8 S6.267 assessment: full nonlinear spatial reduction and retained translations

S6.267 constructs an explicit nonlinear spatial shape and the complete
cotangent momentum lift for the original parent on a local finite-torus
domain. The nonzero spatial constraint complement is solved, while the
three global translation constraints are retained. The same pure nonzero
Gaussian reference mode prescriptions have an exact common zero-translation
charge factor under an explicitly identified finite torus-mode evaluation.
Original P8 remains OPEN.

This is a classical local reduction and a symmetry result for the original
linear-CCR reference, not a nonlinear Gaussian pushforward, support or
interacting evolution theorem. No homogeneous quantum state or continuum
limit is supplied. The full source, both matter spectators, all vector
channels, physical normalization and both historical errata remain.

For Q=det(gamma)^(1/3)gamma^-1, the exact gauge requires div Q=0 and det Q=1.
On the2pi L torus, define P(k)=I-kk^T/k^2 for nonzero k and P(0)=I.
Writing Bf(k)=P(k)f(k), its trace is(2I+Pi0)f: the generated mean has
eigenvalue3, not2. The homogeneous tracefree shape sector has five
components, not two TT polarizations.

In weighted Fourier l1 with weight(1+|k|)^8 and matrix operator norm,
the complete scalar equation for Q=I+tau+Bf is a contraction on
||tau||<=1/100 and ||f||<=1/1000. Its full self-map bound is
366993/10^9 and its Lipschitz constant is67089/10^6<1/10.
These estimates use all convolution products; they are not a truncated
mode determinant. Reality, positivity, exact unit determinant and
zero divergence follow for the entire branch.

The scalar derivative K_Q h=cof(Q):Bh has inverse norm below5/9.
Its inverse and all cofactor contacts enter the full shape tangent.
The exact coefficients f2=(1/2)(2I+Pi0)^-1 trace(tau^2) and
f3=(2I+Pi0)^-1[trace(tau Bf2)-det tau] retain the generated homogeneous
piece. Their full cubic approximation error is at most60|epsilon|^4
for unit-norm direction tau and |epsilon|<=1/100.
Derivative factors2! and3! are distinguished from Taylor coefficients.

The metric is gamma=a^2 exp(2v)Q^-1. At the same linear reference
tau=-t_TT, with the dual sign also reversed. This does not erase
S258's scalar-volume and TT contacts when moving a generic off-slice
metric into this slice. The conformal gauge is exactly the same for
the correct physical metric C gamma even at spatially varying lapse.

The source packet differentiates the complete metric/vector/matter
Lie pairing and retains its full boundary flux. Only the lapse and
normal-vector primary momenta are removed. In particular, the
-W_i partial_j Pi_W^j vector Gauss term remains.

Let A be the metric gauge derivative, V the full Lie action including
both spectators and W, and M_Q=A V the actual nonsymmetric S258
ghost block. On the exact slice its inherited mean-zero weak
coercivity is3/4, with L2 inverse bound4L^2/3.
For every ambient tangent the full horizontal projection is
I-V M_Q^-1 A, not a metric-only correction.

A reduced covector has the explicit ambient extension
pi0=(Pi_v/6)gamma^-1+DQ^*(P_TT Pi_tau), with all matter/vector
momenta retained. If D0=V^*alpha0 is its entire original generator,
the unique lift solving all mean-zero momentum constraints is

alpha=alpha0-A^*(M_Q^*)^-1 Pmean0 D0.

This full correction vanishes on slice tangents, so the complete
canonical one-form pulls back to the physical scalar, shape,
vector and two-matter cotangent form. The direct shape momentum
includes both the K_Q inverse and its adjoint contact; a flat TT
projection of pi would omit them.

The argument uses bounded maps L2 to H^-1 to H1 and back to L2
at fixed smooth coefficients, plus the original derivative-gap
local slice. It does not assume a differentiable same-regularity
diffeomorphism group or an arbitrary strong infinite-dimensional
Darboux theorem. The complete second-class inverse keeps the
potentially nonzero constraint-constraint block and nonsymmetric M.

Constants remain in the ghost kernel. Mean-zero descriptors are
not a Lie algebra: their sin/cos example generates a nonzero
translation. The three residual charges are therefore explicitly

J_i=integral Pi_A partial_i q_A,

summed over all eight nonzero-mode physical configuration channels.
The exact source generator gives this expression with the correct
real Fourier normalization. These three charges commute, and the
cotangent correction leaves them unchanged. A global regular
diffeomorphism quotient or uniform infinite-volume inverse is not
claimed.

The original scalar/tensor Gramian, both symmetric canonical
boundaries, kappa/volume normalization, central time-dependent
swap, exact H SLE and full three-mode Proca graph are all retained.
For a finite real mode pair, a translation-invariant pure Gaussian
has general complex graph Z=[A,B;-B,A], with A symmetric and B
antisymmetric. Odd internal cross blocks are allowed. Its Schwartz
wavefunction is annihilated by the complete configuration-rotation
generator because(Sq)^T Zq=0. The finite compact translation
representation is periodic and its three self-adjoint generators
strongly commute; this is not a generic oscillator phase-space
rotation with a spurious metaplectic sign.

Purity matters: an invariant equal mixture of opposite unit-charge
states has zero charge mean but variance1. Mere invariance of a
mixed density matrix would not show zero-charge support.
No projection, re-normalization or reference-state reset is used.
The finite nonzero-mode factor is kept separate from any unspecified
homogeneous quantum factor and from a nonlinear state on the local
classical chart.

Combining this spatial construction with [S266's full auxiliary
branch](assessment-2026-09-14-p8-affine-nonlinear-lapse-branch.md)
requires that the complete reconstructed twelve bounce invariants
lie inside its actual box. No evaluated canonical-radius or
Gaussian probability for that intersection is asserted.
The full nonlinear and original linear-Gaussian lapse are not
identified.

Independent finite Fourier diagnostics use8K points per coordinate,
so every degree-three product of K-supported fields is non-aliased.
They explicitly measure the full discarded-mode residual rather
than declaring a projected fixed point to solve det Q=1.
Generated means are retained at all cutoffs. The numerical examples
have the stated unweighted A0 norm; they do not claim the same
radius in the stronger A8 norm. The full proof itself is in A8 and
also works in A0.

Further diagnostics test complete metric derivatives, random
nonsymmetric horizontal/dual lifts, nonzero constraint brackets,
all three non-collinear translation charges, general correlated
and squeezed pure covariances, and the mixed-state counterexample.
Nonlinear algebraic fixtures outside the small theorem ball do not
enlarge the proved domain.

The first assembled symbolic run stopped on an unexpanded
real-Fourier normalization residual. Expanding that exact integral
gave zero without changing its formula. The corrected complete
packet run passed in292.136435583001 seconds. All209 independent
numerical/contract tests passed in2.00 seconds with four complete
proof packets deselected.

The complete contract has78 named residuals,1266 scalar entries,
38 proof gates, nine controls and150 rejected inputs, retaining
nine primitive and123 matching rows. Matrix components are not
1266 independent theorems. The written Banach, elliptic, cotangent
and operator-domain arguments are not FORMALIZED.

The final private original-SymPy preflight passed in291.51018845802173
seconds, checking all18 source files,20 AST-derived report fields,
all exact residuals, gates, controls, rejected inputs and frontier
rows. All213 science tests then passed in0.46 seconds with the proof
packets cached. RuntimeWarning was an error. Ruff fixed eight
lint issues, formatted eight files, and its final check and format
check passed. All seven sibling-module export contracts resolved.

All18 scientific files,76237 characters, were frozen at
2026-09-14 16:41:40 UTC only after every source hash matched that
successful preflight. Fresh repository preflight passed in
289.5476244590245 seconds, followed by all213 science tests in
0.58 seconds with the proof packets cached. Its zero exit was
consumed at2026-09-14 16:47:00 UTC. Every private/preflight/repository
source hash matched.

The original-SymPy native rebuild protected5457 scientific inputs before and after construction. Its198813-byte ASCII report was transferred in17 exact chunks with16 observed CONTINUE acknowledgments, ending at offset192000 with6813 characters. The raw un-reserialized SHA256 is36be009c0d58f8e99a4b703bf5a014dde6a28c07dbb888453a8d04aa0fba7f98. Raw validation passed at2026-09-14 16:50:25 UTC for all18 native/private/repository source hashes, all20 AST-derived fields, every count and the complete retained frontier. No chunk was truncated.

The ordinary original-SymPy package replay completed with exit status zero: all238 tests passed in2934.80 seconds (48 minutes54 seconds). Both captured test files were present and unchanged. Its zero exit was consumed at2026-09-14 17:44:13 UTC.

The separate original-SymPy CLI replay completed with exit status zero, consumed at2026-09-14 17:50:43 UTC.

The separately captured full P8 regression completed with exit status zero: all62921 tests passed in5382.45 seconds (1 hour29 minutes42 seconds). It captured787 test files with SHA2560b2bf16f840ca7d77707f85161a78dcfe0dd9c32603c8649b38092129a3c7df0, and all787 were present and unchanged at collection. Its static namespace contained609 ancestors. The independently audited exact-GCD runner passed128 original-path comparisons; final counters were68046 domain fallbacks,7388 exact descents and94 mixed fallbacks. The completed zero exit was consumed at2026-09-14 18:20:54 UTC. The retained full output contains the snapshot, collection, final summary and counters without a truncation marker. The captured snapshot includes S266 and S267, but not the later S268 package.

Native/direct/ordinary/CLI use original SymPy. Only the captured full
regression uses the independently audited exact-GCD adapter and its
original-path selfchecks. Publication requires observed zero exits
for the exact snapshots, native raw SHA/source equality, exact
scoped staging and independent remote verification. No unrelated
P4/P9 changes are included.

The next quantum gate is an explicitly named finite nonlinear
quantization with its entire source, moving-chart time connection,
reference state, ordering and regulator differences retained.
A bounded classical symbol does not by itself imply a positive
Weyl operator, and a Gaussian reference does not have compact
phase-space support. Neither a formal Ward identity nor a
translation-invariant finite diagnostic closes those missing
physical estimates. Original V/G/B/P8 remain OPEN, the completed
scoped P8(a) objective is unchanged, and S261's physical-parent
binding remains explicitly refuted.
