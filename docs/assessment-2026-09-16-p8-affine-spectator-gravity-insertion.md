# P8 S6.305: spectator reference, signed rate and Gaussian response bound

Native, original-SymPy ordinary/CLI and complete captured-suite acceptance
have passed. Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

## Scientific result and prior ownership

S289 already derived the complete crossed source-fixed H/Proca amplitude,
both optical cuts and the local polynomial. S305 independently reconstructs
those results; it is not their first derivation. Its new results are the
uniform full-Born-normalized signed reference rate including massless M1,
the forward Newton endpoint, the Gaussian-only disk16/Dyson remainder
bound and an independent numerical Feynman-branch calibration.

A post-freeze external audit explicitly compared the entire nonlocal,
local and Newton parts with S289 under a common integration variable.
Private prior289-crosscheck.py, session 38771, passed in 5.889832s.
No frozen scientific source or report was edited to add this clarification.

On mu=1,25/4<=s<=16 and every physical nonforward angle, the stated
source-fixed known spectator reference satisfies
-3e-602 < 2Re(M_HV+M_M1)/(A_m+A_G) < 0.
The proof uses the full positive Born denominator, complete crossed
kernels, the original M1 logarithm and the fixed Gaussian Newton shift.
The forward normalized endpoint is -2*delta_kappa_fixed/kappa.

The sign is only that of this known reference. Unassigned physical
finite matching can change the full rate and adds its own Newton
endpoint. It is not a pole-subtracted positivity theorem.
The massless transfer logarithm is retained, not dropped from the
complex amplitude or the Regge problem.

For the same fixed H/Proca Gaussian kernels, the response quotients
have no additional zero in the complex disk |p|<=16. Repeated insertion
of those kernels has a geometric remainder below 5e-1204.
This is not a physical Wilsonian cutoff or a bound on independent
interacting higher-loop graphs.

Independent scalar and full nine-polarization Proca cuts, component
stress contractions, spectral/radial quadratures and complex lower-bank
boundary values calibrate the normalization. The branch-safe closed
radial function uses atanh(sqrt(p/B))/(B*sqrt(p/B)); replacing its
denominator by an independently chosen principal sqrt(p*B) is rejected.

Physical hard matching, interacting quantum state, all-loop/all-N
errors, finite forward cross section, complex Regge control, the
common-parent bounce and full V/G/B/P8 closure remain outside the claim.

## Acceptance evidence

Cold 48147 passed the exact/export/source/serialization preflight in
86.300393s and all 395 scientific tests in 0.90s.
Independent fresh 96319 forced the entire immediate S304 parent's
packets, passed in 88.445194s and passed all 395 scientific tests in 0.97s.
All three copies of the 18 ASCII sources and their hashes agree.

The frozen sources total 67,343 bytes, with 218 named/scalar checks,
44 proof gates, 8 controls, 73 rejected inputs, 9 primitive records,
161 matching records and 6 historical qualifications.
Native helper 75613 received one specification and protected 6,164 inputs.
Its 183,181-byte ASCII report arrived in 16 chunks with 15 acknowledgments.
SHA256: 22b7a327ff4d4a2c39fc6ef4fe7435e84549def254317347e7f8abc01a362011.
Both raw/source/count/frontier validators and canonical 20-field REPORT
AST passed. Raw freeze: 2026-09-16 03:49:55 UTC.

Original-SymPy ordinary 65879 passed all 420 tests in 3283.11s (0:54:43).
Independent CLI 59201 returned C0. Both started
2026-09-16 14:49:52 UTC; completion was captured
2026-09-16 15:44:49 UTC. Complete outputs are retained.

Shared FULL 78849 passed 80,430 tests in 5795.81s (1:36:35), with all
861 captured files unchanged. Snapshot SHA256:
a5d9f80fc3d21b57f0728c407d294b382cba7a927d24196171d12a442b7e7329.
It captured 683 namespace ancestors, 5 helpers and 128 adapter contracts.
Counters: 68,312 domain fallbacks, 7,388 exact descents, 94 mixed fallbacks.
It includes S304/S305, not S306 or later. Completion was captured
2026-09-16 05:38:42 UTC; the complete output is retained.
The later 81,447-test FULL 71089 also passed through S307.

Native/direct/ordinary/CLI retain original SymPy; the audited exact-GCD
adapter is FULL-only. Publication rechecks every frozen hash, raw report
and the exact 22-file scope, excluding unrelated P4/P9 work.
Original V/G/B/P8 remain OPEN.
