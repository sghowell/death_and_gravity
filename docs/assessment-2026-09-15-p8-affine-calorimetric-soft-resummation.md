# P8 S6.299: complete leading-soft calorimetric resummation

The whole leading-soft total-energy sum is explicit, including its finite
gamma-function factor and noncommuting regulator/detector limits.
Native replay, original-SymPy ordinary/CLI, and the captured complete P8
regression passed. Original V/G/B/P8 remain OPEN.

## Scientific result and boundary

For the inherited physical massive-scalar soft index a and finite
real/virtual conversion Delta, the complete leading-soft factor is

exp(Delta) exp(-gamma_E*a) x^a / Gamma(1+a).

The total emitted energy is constrained by one calorimetric cutoff.
Its N-particle simplex and the entire Poisson sum are retained before
the infrared regulator is removed. A fixed-N limit is not substituted
for the complete sum. For0<x<=1 the physical positive Poisson cumulative
distribution is exp(-gamma_E*a)x^a/Gamma(1+a); exp(Delta) converts to the
stated inherited analytic virtual reference, not a different probability law.

At the unchanged original parameters the finite prefactor differs from
1+Delta by less than20000/kappa^2=2*10^-1596. The power x^a is unexpanded:
there is no unjustified uniform Taylor approximation when resolution
becomes exponentially small. The x=exp(-kappa*chi) detector scaling has
its explicit limiting factor, while a simultaneous regulator scaling
can give a different limit. The admissible order of limits is recorded.

This is the complete leading-soft model, not a bound on arbitrary-N
nonleading radiation, hard-loop corrections or a quantum UV completion.
Radiative-state soft-index changes, hard evanescent terms, matching,
complex Regge and the original state/bounce obligations are not set to
zero. Scoped P8(a) and all frozen ancestor qualifications are unchanged.

## Acceptance evidence

Cold preflight24319 passed in89.881584833s and all373 science tests
passed in0.65s. Fresh-copy67168 forced the whole immediate parent
packets and passed in92.364864500s;373 science tests passed in0.67s.

The18 ASCII sources total55,714 bytes:223 named identities,223 scalar
entries,37 proof gates,8 controls,73 rejected inputs,9 primitive records,
155 matching records and6 unchanged historical qualifications.
Private/fresh/repository source bytes and hashes agree.
Source freeze:2026-09-15 23:38:02 UTC.

Native helper75613 received the specification once and protected6,050
inputs. The154,347-byte ASCII report arrived in13 chunks with12
acknowledgments. SHA256:
604445d8362c3dd9cec6e456bab0743e0611aaac2e60794e64f99be41cd63879.
Both raw/source/count/frontier validators passed, including the20-field
and canonical REPORT-path AST checks.
Raw freeze:2026-09-15 23:42:42 UTC.

Original-SymPy ordinary45152 passed all398 tests in 3237.44s
(0:53:57), with both captured files unchanged. Independent CLI57295
passed. Both started2026-09-16 02:06:47 UTC; completed acceptance was
captured 2026-09-16 03:02:45 UTC.

Shared FULL89263 passed 78,218 tests in 5702.41s (1:35:02), with all
851 captured files unchanged. Snapshot SHA256:
b6f53a55404237f6880731ab58e77b368a7718c1be1c81626b6d9f0c6fbb0444.
It captured673 namespace ancestors,5 helpers and128 adapter contracts.
Counters:68,283 domain fallbacks,7,388 exact descents,94 mixed fallbacks.
It includes S298/S299/S300, not S301 or later.
Completion was captured2026-09-16 02:07:49 UTC; complete output retained.

Native/direct/ordinary/CLI retain original SymPy. The exact-GCD adapter
is FULL-only. Publication rechecks all frozen hashes and the exact22-file
scope, excluding unrelated P4/P9 work.
Original V/G/B/P8 remain OPEN.
