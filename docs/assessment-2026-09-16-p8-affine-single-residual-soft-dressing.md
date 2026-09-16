# P8 S6.309: state-correct soft dressing of one finite residual

Native, original-SymPy ordinary/CLI and complete captured-suite acceptance
have passed. Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

## Result and exact scope

The unchanged source has mu=nu=1, kappa=1e800, 25/4<=s<=16,
every physical nonforward hard angle and 0<x<=1/8. The normalization
is the positive complete Born amplitude, including matter and gravity.

The seed is the finite signed one-real residual from all 47 radiation
trees: 26 matter and 21 Einstein trees, with all their interference.
For its cumulative total variation, S307 gives
V_R(x)<min(1e32,2e14*x+2e37*x^2)/kappa.

This checkpoint dresses exactly one such residual with arbitrarily many
additional LEADING-soft gravitons. The complete massive-plus-marked-null
current of the same radiative state determines both its real soft index
and its virtual pole. The detector cut leaves y=x-omega for the additional
emissions; the Bose relabeling gives exactly 1/N!, not an extra marking
factor. The resulting named reference is

D[R](x)=integral_(0<omega<x) P_sigma(x-omega) dR(sigma),
P_sigma(y)=exp(Delta_sigma)*exp(-gamma_E*a_sigma)
           *y^a_sigma/Gamma(1+a_sigma).

The fractional-dimensional angular contraction can be negative at
finite regulator. A uniform absolute bound on the entire signed series,
including the zero-physical-index boundary, justifies the limit against
the signed seed. No fictitious finite-regulator probability law is used.
The regulator is removed at each fixed positive detector threshold;
no arbitrary joint regulator/forward/resolution limit is asserted.

Relative to the unexpanded elastic soft reference P0,
|D[R](x)|/P0(x)
 <(1+4250/kappa)*min(1e32,2e14*x+2e37*x^2)/kappa
 <2e-768,
and separately
|D[R](x)|/P0(x)<(4e14*x+4e37*x^2)/kappa ->0.
Both bounds are uniform in every nonforward hard angle. Detector-energy
powers and the remaining-energy endpoint are retained exactly.

Thus P0+D[R] is positive pointwise on the stated domain. This does NOT
establish monotonicity, a positive detector measure, unitarity, or equality
to the complete interacting rate. Additional emissions in this reference
are leading soft with the marked state fixed. Multiple finite residuals,
simultaneous soft contact terms, correlated exact recoil, finite radiative
hard loops, evanescent hard terms and physical matching need separate
estimates. No omitted coefficient is selected or set to zero.

## Independent evidence

Independent diagnostics check marked Bose counting through four labels,
literal one- and two-soft simplices, marked energy moments, and the
complete finite-regulator series at decreasing regulator. Positive-index
and zero-index cases include an alternating signed boundary.

Convolution calibrations keep a state-dependent index and finite
prefactor, signed seed, remaining-energy endpoint and an exponentially
small threshold. Three actual original 47-tree samples independently
check negative residual densities, positive full Born and the complete
massive-plus-null index. These calibrations are not certified numerical
quadrature error estimates and do not replace the written uniform proof.

The first scratch convolution input used a negative index slope that
crossed below the stipulated physical-index domain. Changing only that
unfrozen synthetic slope from -10*omega/kappa to -2*omega/kappa restored
the domain. Assertions, tolerances and backend were unchanged; the
corrected complete probe, including actual-source samples, passed.

Cold 45641 passed the complete original-SymPy preflight in 88.837986s
and all 577 science tests in 2.75s. Independent fresh 71722 forced the
entire S308 parent packets, passed in 92.101490s, and then passed
577 science tests in 2.65s. All 18 private/fresh/repository sources were
byte-identical before freeze.

The 18 ASCII sources total 70,532 bytes: 352 named/scalar checks,
63 proof gates, 8 controls, 73 rejected inputs, 9 primitive records,
165 matching records and 6 historical qualifications.
The native helper received one specification, protecting 6,240 inputs.
Its 177,949-byte ASCII report arrived in 15 chunks with 14 ACKs.
SHA256:
22e89710b44e3df51ab35cf4ef541f70c42c68ee740a29f261f9f9d16694b59c.
Both raw/source/count/frontier validators and the 20-field REPORT AST
passed. Native raw freeze: 2026-09-16 06:45:05 UTC.

## External acceptance

Original-SymPy ordinary 88284 passed all 602 tests in 3258.84s (0:54:18).
Independent CLI 18228 returned C0 with the scoped success verdict.
Both lanes began approximately 2026-09-16 18:44:20 UTC; their successful
completions were captured 2026-09-16 19:39:03 UTC.

Shared FULL 40862 passed 82,577 tests in 5644.89s (1:34:04), with all
869 captured files unchanged. Snapshot SHA256:
f0fcbc72c5385dca25337e1a08dbac522fc16f197223daa3208119cb55f18bfa.
Its 691 namespace ancestors, 5 helpers and 128 adapter contracts were
captured; counters were 68,421 domain fallbacks, 7,388 exact descents
and 94 mixed fallbacks. It includes S308/S309, not S310 or later.
Completion was captured 2026-09-16 16:25:14 UTC.
The later FULL 50300 also passed through S311: 83,650 tests in 5799.26s,
873 unchanged files, snapshot SHA256
bad08d33f93588b01f48977c40f9fd9ea2b9f8345f6966f1a4817653528931b9.

Native/direct/ordinary/CLI retain original SymPy; the audited exact-GCD
adapter is FULL-only. Publication checks the exact 22-file scope,
frozen source/raw hashes and preservation of unrelated P4/P9 work.
Original V/G/B/P8 remain OPEN.
