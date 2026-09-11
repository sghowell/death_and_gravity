# P8 continuation: conditional affine-Proca state and clock stress

S6.176 keeps the same CD-REG-AFFINE-ISO classical parent.
Original P8(b), V/G/B and P8 remain OPEN.

## What is established

Specialize the allowed vector parameter to zeta=1e-6, with kappa=1e800.
The EXACT canonical field and source are

    A_mu=sqrt(kappa*zeta)W_mu,  J_mu=sqrt(kappa/zeta)S_mu,

so the conditional Gaussian sector is ordinary mass1000 Proca with its
complete source and source-squared contact. W does not have the
unscaled canonical covariance. Only the three physical vector modes
are quantized with the explicitly stated canonical measure; the
full affine/metric/light/auxiliary measure is not supplied.

The arbitrary-lapse and scale-factor action is varied before its
temporal constraint and canonical substitution. The physical energy
and pressure, longitudinal normalization and full clock frequencies
therefore match the source-pinned ordinary-Proca calculations as
operators and readouts, not merely as oscillator equations.

The same all-order cutoff-summed Cauchy covariance from S6.55 is used
at t0=-1/2; the older finite-order WKB reference is not called Hadamard.
For compact smooth histories after a common Cauchy neighborhood,
the sourced mean is the retarded Proca solution and the connected
Hadamard covariance is unchanged. The coherent field is a state of
the SOURCED affine algebra, not a false homogeneous bisolution.

The Proca time-slice/Hadamard propagation hypotheses are checked.
The original causal-wavefront gap in the Moretti/Murro/Volpe argument
is handled using Fewster's corrected Theorem5.1, not ignored.
These external results supply the microlocal state property; they
do not supply any numerical stress constant.

## Covariant prescription and quantitative result

The specified dimensional local prescription is the constant-mass
ordinary-Proca specialization. Its nonminimal normal-clock terms
vanish identically as functions, including off the clock, so the
remaining covariant curvature densities have no hidden X denominator.
The finite heat-coefficient terms, including the dimensional scalar
rank contribution, are kept. No arbitrary normal-ordering or old
scalar stress-canceling profile is installed.

On the unchanged clock and I=[-1/2,1/2], the exact full-continuum
state, subtraction and remainder bounds give all twelve estimates

    |partial_t^j rho|, |partial_t^j P| < 1e30,  0<=j<=5,
    |partial_t^j rho|/kappa, |partial_t^j P|/kappa < 1e-770.

These are absolute specified conditional-vector contributions.
The normalization is the fixed kappa, not the zero classical density
at the bounce. A small prescribed-clock residual is not a solved
quantum-corrected background. The numerical constants are confined
to the stated slab, not asserted uniformly over all real time.

The full clock source and its first variation vanish, while its
second variation can be spatially nonclosed. The sourced mean
starts at second order, its free classical stress at fourth order,
and the full source/contact force at third order. Temporal
constraints, seagulls and higher source insertions cannot be dropped.
C5 time-derivative bounds are not functional metric-response or
noise norms.

## Verification

Final private preflight passed, followed by207 science tests in2.70
seconds with the computed packets cached. Fresh repository science
passed207 in156.40 seconds. Ordinary replay passed232 in2119.04
seconds, and the independent native CLI replay passed.

The complete captured P8 snapshot passed30234 tests in3335.66 seconds,
final exit code0. All605 captured test files remained unchanged.
Path-list SHA-256:

    d1cf35580f712b58631909c183e3abee0c19d8890ec50ede3945b13945fdbe58

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34177 domain fallbacks,7340 exact descents,88 mixed fallbacks.
Native, direct science, ordinary and CLI used original SymPy.

The native74681-character report was transferred losslessly in seven
chunks and independently checked against all18 frozen scientific hashes.
Its20 fields contain64 named identities,94 scalar entries,19 proof
gates,9 controls and100 rejected inputs. Report SHA-256:

    1175702d6980db090396863b807af0e767241daeb7d8dabf138ed11c67e1ae75

The report fully rebuilds S6.175 and additionally checks/rebuilds the
explicit S6.55 state and S6.82 stress sources. Independent controls
vary full lapse/scale actions, reconstruct constrained frequencies,
check actual Cauchy preparation, physical curvature-density variations,
the local Ward identity and nonlinear source contacts. A private
wording-control mismatch was corrected before freezing.
These are exact algebra and written proofs, not FORMALIZED.

No frozen byte was changed. The exact22-file checkpoint excludes
S6.177, the private S6.178 draft and unrelated P4/P9 work.

## Next frontier

S6.177 has passed repository/native checks and is in fresh replay.
It constructs an actual gravitational-decoupling family through
the complete base action, preserving all independent canonical
scalar functions and vacuum masses. The finite-gravity forward
pole is kept separate from a fixed-transfer limit.

S6.178 is a private quantitative retarded-response draft. Its initial
science tests pass for a finite-time source-energy estimate, the full
near-clock source and a cubic fixed-metric light-force bound, with
homogeneous cancellation and spatially nonclosed controls.
It is not yet a published or fully replayed certificate.

Full metric/noise response, interacting measure and state, light/
tensor/auxiliary/mixed loops, physical cutoff, controlled background,
vacuum contour/cuts/truncation and finite-gravity IR/Regge estimates
remain open. Scoped P8(a) and A.20-A.23 are unchanged.
No user-intervention blocker has been identified.
