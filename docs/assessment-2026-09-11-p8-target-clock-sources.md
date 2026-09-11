# P8 continuation: actual target clock sources and constraints

S6.173 identifies the leading terms a common parent must reproduce.
Original P8(b), V/G/B and P8 remain OPEN.

## Full target, not a vacuum-only coefficient match

The unchanged S6.109 analytic family has the original CD_matter clock
jets, with the retuned margin and original nonzero M1 momentum kept.
The actual geometry is a=(1+t^2)^2. The physical action is kappa times
the dimensionless kernel, with canonical clock Phi=sqrt(kappa)u.

The covariant action is varied with arbitrary lapse before setting
N=1. The nonminimal curvature coefficient vanishes in value on the
clock but its X derivative contributes. Likewise the A3 action value
vanishes there while its lapse-velocity variation contributes.

At the bounce the null Euler contributions, divided by kappa, are

    Einstein:                  +8
    canonical clock kinetic:   +1
    F minus that kinetic:      -2101/100
    nonminimal curvature:      +24
    A3:                        -12
    original M1:               +1/100.

Their sum is exactly zero. The full energy and pressure equations
cancel at every real time, with original matter conservation retained.
Individual action pieces are not separately conserved matter tensors.
The canonical reference here is the kinetic X/2 only; the entire lower
potential remains in F minus that reference.

A4 and A5 have zero first background variation, but they cannot be
removed from the action: they supply the trace/lapse velocity-Hessian
degeneracy. Deleting them changes the rank and produces determinant
-(4h+3)/(3h^2), rather than the correct rank-one perfect square.
This is a structural control, not a new propagating ghost classification.

The matter field must also be retained during metric variation.
Substituting its conserved-momentum solution into the Lagrangian too
early reverses the energy/pressure signs. The fixed-momentum Routhian
restores the required Legendre term and correct signs.

The old vacuum action estimate uses the fixed Fourier-unit real
Schwartz class, whose first derivatives have supremum at most1.
The canonical clock derivative is1e400. Its normalized first-jet
distance from that class is at least1-1e-400, even on a finite interval.
A small action norm on the vacuum class therefore does not control
clock equations, their variations or a common parent solution.

## Verification

Final private science passed165 tests in4.03 seconds.
Fresh repository science passed165 in4.14 seconds.
Ordinary replay passed190 in2019.11 seconds; independent native CLI
replay passed. The complete captured P8 snapshot passed29526 tests
in3333.98 seconds, final exit code0.

All599 captured test files remained present and unchanged.
Path-list SHA-256:

    75bf319e6e69cfbbe1a392e7660eab06dedfd3f0ef7599f917ba3b5bd0f49188

The full-only exact GCD adapter passed128 original tuple comparisons.
Final counters:34090 domain fallbacks,7340 exact descents,88 mixed
fallbacks. Native, direct science, ordinary and CLI used original SymPy.

The native28186-character report was transferred losslessly in three
chunks and independently checked against15 scientific hashes.
Its20 fields record57 named identities,57 scalar entries,18 proof gates,
9 controls and55 rejected inputs. Report SHA-256:

    3740dfeae08879dda0039c10fd94c09018ef964d1ca12690624089927a852ba5

Independent controls include direct connection/scalar contractions,
source-pinned clock jets,14 compact lapse/scale functional variations,
source omissions, velocity ranks, matter Routhian variations and
normalized matching-domain distances. The written analytic statements
are not claimed as FORMALIZED or as a quantum background theorem.

No frozen source or report byte was changed. Exact19-file staging
excludes later continuations and unrelated P4/P9 work.

## Next frontier

S6.174 now provides a separately named classical CD-REG-AFFINE-ISO
candidate on the same analytic domain. Its native report, all source
hashes and repository science checks passed; fresh ordinary, CLI and
full replays are running. It uses a regular metric chart, a new
source-centered isotropic trace mass and a vacuum-regular local shift.
Its source vanishes through first order on the exact clock, and the
joint lapse/vector Jacobian is nonzero at every finite clock time.

This is a classical common-domain construction, not a transferred
GY14 quantum parent or a completed B gate. Its nonlinear vector
source, real-time response and quantum stress require explicit
treatment. Physical cutoff, omitted-order bounds, quantum matching,
controlled background/response, vacuum contour/truncation and
finite-gravity IR/Regge estimates remain open. Scoped P8(a) and
A.20-A.23 remain unchanged. No user-choice blocker has been identified.
