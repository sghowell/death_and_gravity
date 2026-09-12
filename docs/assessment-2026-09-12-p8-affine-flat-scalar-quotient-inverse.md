# P8 continuation: full flat scalar quotient reference inverse

S6.224 completes a two-channel scalar gauge-quotient reference inverse, uniformly in all spatial momenta and without spatial derivative loss. It retains the [original isolated shear pole](assessment-2026-09-12-p8-affine-isolated-shear-resolvent.md). It is not yet the full curved quantum inverse or a solution of the [S222 force graph](assessment-2026-09-12-p8-affine-quantum-forced-constraints.md). Original P8 remains OPEN.

## Full covariant projectors and both ordered gauge legs

For q=P^2, p=lambda^2+q and spatial scalar amplitudeQ=2zeta I-2c Pi, define

S=(D^2+2q/3)zeta-D^2c/3,
W=-q zeta-D^2c,
B=[[D^2+2q/3,-D^2/3],[-q,-D^2]].

The full S200 spatial covariant projectors give p^2 Pi0(QD,QG)=12SD SG and p^2 Pi2(QD,QG)=(8/3)WD WG. These are not center-of-mass projectors applied at nonzero transfer. The original finite Hessian supplies-4SD SG-(4/45)WD WG. The complete normalized quotient reference is

A_q=B^T diag(Ftrace(D^2+q),(8/3)F2(D^2+q))B.

The three-source reference in(n,zeta,b) has a gauge kernel. Its curvature map is

M(lambda)=[[q/3,lambda^2+2q/3,-lambda/3],[q,-q,-lambda]].

The source gauge vector is(lambda,0,q); the detector has the opposite time sign. Writing n=D eta, b=Dc+q eta, the complete source map givesM T=[0,B], and the separately transformed detector leg gives the corresponding adjoint reduction. The original three-source operator isM(-lambda)^T diagF M(lambda), not a naively symmetric fixed-lambda matrix. No inverse of that gauge-degenerate three-source block is claimed.

## No inverse-transfer singularity in the causal coordinates

The exact B inverse is

[[1/p,-1/(3p)],
 [1/p-1/lambda^2,-1/(3p)-2/(3lambda^2)]].

The identity q I^2 L_q^-1=I^2-L_q^-1 avoids division byq. Its causal kernel is

R_q=[[s_q,-s_q/3],[s_q-t,-s_q/3-2t/3]],
s_q=sin(sqrt(q)t)/sqrt(q), s_0=t.

R_q(0)=0 and ||R_q'||<5/2 uniformly inq. Indeed its exact Frobenius square is(20cos^2-14cos+13)/9<=47/9<25/4. Thus ||R_q(t)||<5t/2. The continuousq0 amplitude extension is not a substitute for solving literal homogeneous lapse/shift constraints.

## Both channels, the retained pole and an ordinary composite inverse

The shifted trace and shear inverse spectral measures use frequencies sqrt(q+tau), with the shear pole atsqrt(q+r m^2) retained. Their positive static masses are bounded by1/4 and30, so the primitives satisfy

||Jtrace,q||<=1/2, ||J2,q||<=60,
||Jdiag,q||<=45/2, Jdiag=diag(Jtrace,(3/8)J2).

The complete inverse is the ordinary matrix convolution

E_q=R_q' * Jdiag,q * R_q^T.

Its Laplace transform isB^-1 diag(1/Ftrace,(3/8)/F2)(B^T)^-1, and both products withA_q are the identity. Causal distribution definitions retain every initial boundary term; no initial atom or shear pole is removed. The source-pinned forward factors are defined by paired local/cut distributions, avoiding divergent separate pieces.

The simplex bound gives

||E_q(t)||<=375t^3/16,
||E_q||<=375T^4/64<6T^4

on bothC_tH^r andL2_tH^r, for every realr and ALL comoving momenta. No spatial derivative loss or transfer cutoff is used. The first-time-derivative kernel is bounded by1125t^2/16 and its operator norm by375T^3/16.

The complete causal forward graph includes the initial boundary. The inverse is bounded on ordinary source spaces, while the forward distribution need only have a finite-order mapping into a weaker space; no bounded forward self-map or full curved graph invariance is asserted.

## Keep physical normalization and scope

The amplitude metric Gram matrix is[[12,-4],[-4,4]], strictly between2I and14I. Physical metric/source norm conversions must retain those constants. The displayedA_q is64pi^2 times the amplitude-coordinate Hessian. Restoring physical Hessian units multiplies the inverse by64pi^2. The flat unit-density S222 force reference also restoreskappa, giving the bound375pi^2 kappa T^4. This is not a kappa-small full feedback estimate.

The actual curved mass/state/contact/tree/matter remainder and its compatible-space bridge are not supplied by a flat reference factorization. Neither quantum stability nor a controlled background, finite-coupling/nonlinear parent remainder, physical cutoff or original V/G/B closure follows.

## Independent checks and completed fresh validation

Three independently rotated full four-dimensional metric/curvature/projector fixtures check the complete normalization. Further tests separately transform both gauge legs, directly solve the rational matrix systems, solve forced time-domain coordinate equations atq0 and nonzero transfers, and reconstruct both shifted continuum dispersions at high precision with the shear pole included.

An initial independent fixture used Python floating division in nominally exact momenta. Coercing its momentum/frequency inputs to exact rationals fixed two checks without changing the formulas. Corrected independent science passed19 tests in4.37 seconds; the integrated suite passed270. Final private science passed270 in4.97 seconds, repository science270 in5.22 seconds. Preflight verified17 inputs and20 fields in1.62 seconds. Formatting and lint passed.

Fresh original-SymPy ordinary replay passed295 tests in2338.75 seconds. Standalone CLI replay passed. Full P8 regression passed45498 tests in4460.77 seconds, final exit code0. Its701-file snapshot SHA is
`3227125e9cced6813ce0c33665a32f8e5fb89fadcbb9b8853916c58b54384f81`.
The full launcher retained the S219 helper-directory PYTHONPATH allowance. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36072 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native/direct/ordinary/CLI retained original SymPy.

The41967-character native report was transferred in four checked chunks. Its20 fields contain44 named identities,105 scalar entries,25 gates,nine controls and180 rejected inputs. Every frozen source/proof/test hash and report SHA were verified:
`a2287a5083d5dceab6886f4f7de3304ce3a7af56c772e76dc869115ad8ac3fb5`.

The exact21-file publication manifest contains17 inputs, report, this audit, CLAIMS and README. The pending curved and weighted checkpoints, private finite-local research and unrelated P4/P9 changes are excluded. The causal/functional arguments are written proofs, not FORMALIZED.

## Next work

The actual curved-coordinate and weighted bounded-correction successors are under fresh validation. The complete local remainder is being investigated through an exact curvature-coordinate factorization. None yet establishes the complete actual curved quantum-force inverse or original P8 closure. All primitive and prior matching statuses remain unchanged, and continuation has no user-intervention blocker.
