# P8 S6.297: first-Newton inclusive assembly and forward-domain boundary

The complete formal first-loop, first-Newton rate contribution is now
assembled in the unchanged matter on-shell prescription and the S296
analytic-soft reference. Native replay, original-SymPy ordinary/CLI
checks, and the complete captured P8 regression passed.
Original V/G/B/P8 remain OPEN; scoped P8(a) is unchanged.

## Scientific result and boundary

The full finite-kappa, no-internal-graviton one-loop matter diagrams
equal the S239 on-shell diagrams. This follows from the current complete
canonical source functions and the one-loop vertex-count identity,
not by dropping numerically small kappa-dependent coefficients.

With positive matter and gravity Born amplitudes A_m and A_G, the
first-loop rate contribution linear in the Newton marker is

N_1 = A_m^2 [H + Delta_soft + R_real] + 2 A_G Re(L_m).

H retains the independent finite gravitational hard matching. The
selected matter-Born bracket is the S296 real-virtual pairing.
The additional matter-loop interference is finite at fixed nonforward
angle. Against the explicitly unexpanded reference (A_m+A_G)^2,
its absolute relative contribution is below 10^-199 for original
parameters, 4<s<=10^196, and every nonforward physical angle.
Combining it with S296 on 25/4<=s<=16 gives a known conversion and
interference error below 2*10^-199. If L_m is retained explicitly,
the conversion/recoil error alone remains below 10^-792.

An original-parameter forward cone makes A_G dominate A_m despite
kappa=10^800: for w=1-z^2<=10^-204 the ratio exceeds 150, while
w>=10^-190 bounds it below 1.7*10^-9. The normalized error estimate
therefore does not imply a uniform expansion about matter Born.

For the forward dispersion question, the finite contour and graviton
pole must be combined before the forward limit. A written sufficient
Cauchy estimate makes the missing absolute remainder bound explicit;
relative Regge dominance alone does not supply that bound. Its
model-specific analyticity, disk, and norm hypotheses are not proved.

This is a formal coefficient and rate statement, not a UV-complete
S-matrix. Pure-gravity loop terms, gravity-Born radiation, further
Newton orders, finite hard matching, high-energy/forward estimates,
and the original state/domain/measure/common-parent bounce remain open.
No finite matching coefficient is chosen and no frozen source changes.

## Acceptance evidence

Cold preflight 6269 passed in 84.614572208s; all 360 science tests passed
in 0.12s. Fresh-copy preflight 41194 forced the full immediate parents
and passed in 85.924959666s; all 360 science tests passed in 0.12s.

There are 18 ASCII scientific sources totaling 63,118 bytes, 196 named
identities and scalar entries, 34 proof gates, 8 controls, 73 rejected
inputs, 9 primitive records, 153 matching records, and 6 unchanged
historical qualifications. Private/fresh/repository bytes and hashes
agree. Source freeze: 2026-09-15 22:28:54 UTC.

Native helper 75613 received the specification once, protecting 6,012
inputs. The 159,240-byte ASCII report arrived in 14 chunks with 13
acknowledgments. SHA256:
6206f5d56959f01a8467e2c2b74ff43a2fad3cc75df0cee67e9bb8ce82f39564.
Both private and repository raw/source/count/frontier checks passed,
including 20-field and canonical REPORT-path AST checks.
Raw freeze: 2026-09-15 22:31:36 UTC.

Original-SymPy ordinary 24752 passed all 385 tests in 3148.16s (0:52:28),
with both captured test files unchanged. Independent CLI 18038 passed.
Both started 2026-09-16 00:17:44 UTC; completed acceptance was captured
2026-09-16 01:10:51 UTC.

Shared FULL 10410 passed 77,019 tests in 5603.61s (1:33:23), with all
845 captured files unchanged. Snapshot SHA256:
540f6ad544eab8bbe0f03674c3bc2f8b174b4752bed8f183f7e67c5510c1d55b.
It captured 667 namespace ancestors, 5 helpers and 128 adapter contracts.
Counters: 68,291 domain fallbacks, 7,388 exact descents, 94 mixed fallbacks.
It includes S296/S297, not S298 or later. Successful completion and
snapshot checks are retained; an intermediate progress-dot segment
was truncated by the terminal transport.

Native, direct, ordinary and CLI use original SymPy. The exact-GCD
adapter is FULL-only. Publication rechecks the exact 22-file scope and
all frozen hashes. Unrelated P4/P9 files are excluded.
Original V/G/B/P8 remain OPEN.
