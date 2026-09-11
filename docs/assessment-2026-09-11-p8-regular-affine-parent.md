# P8 continuation: regular classical affine parent on the analytic target

S6.174 constructs a separately named CD-REG-AFFINE-ISO classical EFT.
Original P8(b), V/G/B and P8 remain OPEN.

## What changed

The complete S6.109 analytic light target and original M1 are retained,
with kappa=1e800 and the unique vacuum-regular lower dictionary q.
This is a new classical affine/vector parent, not a transfer of the
earlier kinetic parent or GY14 quantum matching.

On the full stated coefficient domain, every real u and
-1/4096<X<6/5, the physical metric admits an explicit regular chart:

    g_physical=C g_hat+D du du,
    C=R^(-1/2), D=(1-C)/X, C+DX=1.

The inverse, metric determinant ratio C^3 and ten-dimensional
fixed-du Jacobian C^9 are nonzero. The coefficient quotients are
analytic through X=0, including nonzero null gradients. The physical
matter metric remains g_physical, not the hatted cone.

A new source-centered isotropic trace mass is checked against all64
connection directions and the full60-dimensional projective quotient.
It leaves56 complementary directions exactly algebraic and retains
four vector components. Both the linear source and its constant
centering term are necessary. The new-to-old quotient determinant
ratio is -tau*sigma^3, not its reciprocal; the first private dense
controls caught and corrected that draft error before freezing.

The exact retained extra action is

    kappa integral sqrt(-g) [-zeta F(W)^2/4+(W-S)^2/2],
    W=T-B(u,X)du,  0<zeta<=1/2000.

The full source S is regular. The lower q cancels from it but is not
reset to the old singular clock boundary. In the timelike hatted
chart its normal component is (R-1)(Khat-3H*sqrt(X)).
It vanishes through first order on the full CD clock, so the classical
clock equations are exact and the light quadratic block is unchanged.
The nonlinear source is not discarded: it generally has a nonzero
second clock variation and can be nonzero at a null gradient.

## Constraints and healthy quadratic blocks

All lapse-time and spatial boundary terms are kept. The ten-velocity
block is nonsingular in the stated range. Joint lapse/vector temporal
elimination gives a clock auxiliary Jacobian diag(-2J,-1), with J a
strictly positive rational function whose numerator coefficients are
all positive. The independent frozen light-principal comparison and
retuned margin agree. Dropping the lapse boundary reverses the sign
of the bounce control and is rejected.

The implicit-function argument gives local seven-mode constraint
neighborhoods at every finite clock time. It does not supply a uniform
neighborhood radius at infinite time, nonlinear stability or a
whole-domain characteristic theorem. The constant vacuum has the
healthy two tensor, one clock, one matter and three massive-vector
quadratic modes. On the clock the normalized vector frequencies have
an explicit floor at least q+1985. This is not a physical cutoff or a
stationary scattering gap.

## Verification

Final private science passed193 tests in8.02 seconds.
Fresh repository science passed193 in8.38 seconds.
Ordinary replay passed218 in1974.12 seconds; independent native CLI
replay passed. The complete captured P8 snapshot passed29744 tests
in3294.26 seconds, final exit code0.

All601 captured test files remained present and unchanged.
Path-list SHA-256:

    76408e101689b9a53cabda30a0c00c7d3f5e53aaf6b0af3c7f570b127c69598d

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34135 domain fallbacks,7340 exact descents,88 mixed fallbacks.
Native, direct science, ordinary and CLI used original SymPy.

The native26073-character report was transferred losslessly in three
chunks and checked against all18 scientific hashes. Its20 fields
record67 named identities,844 scalar entries,19 proof gates,9 controls
and84 rejected inputs. Report SHA-256:

    e0cf546e77cf12ff4e0f2e6bef706ff6de1eddc40cacff3b0b151fcb7009a4f4

Independent controls include dense quotient solves, literal metric
Jacobians, full-function metric variations, null-gradient limits,
vacuum jets, constraint inertia and all-momentum normalization formulas.
The claims are exact algebra plus written proofs, not FORMALIZED.

No frozen source or report byte was changed. Exact22-file staging
excludes later continuations and unrelated P4/P9 work.

## Next frontier

S6.175 has passed native and repository science checks, with fresh
ordinary, CLI and full replays running. It retains the complete source,
derives closed-source/homogeneous classical lifts and the exact Gaussian
Schur identity, and separates a complete-source Euclidean bound from a
leading-germ Lorentzian estimate. Full real-time source control is not
inferred across the physical vector pole.

The next calculation checks a conditional sourced Proca Hadamard state
family and the exact operator, measure and stress bridge at the allowed
m=1000 specialization. This does not transfer an interacting parent,
countertune the classical target or assert a quantum-corrected bounce.

Physical cutoff, omitted-order bounds, quantum matching, controlled
full background/response, vacuum contour/truncation and finite-gravity
IR/Regge estimates remain open. Scoped P8(a) and A.20-A.23 are unchanged.
No user-intervention blocker has been identified.
