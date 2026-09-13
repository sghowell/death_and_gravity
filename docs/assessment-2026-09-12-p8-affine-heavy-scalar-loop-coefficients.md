# P8 continuation: complete first-loop low-energy coefficient matching

S6.237 bounds the three low-energy coefficients b20,b21,b40 of the unchanged separate V2S-T1-OS4 amplitude. It uses a directly proved joint complex neighborhood, not an inference from the real physical-angle estimate. No new derivative counterterm or heavy-weight subtraction is introduced. Original P8 remains OPEN.

## Complete logarithmic reduction through four inverse-mass orders

Keep every ordered bubble, triangle and box in S235 and the complete primitive T from S236. For w=1/n and H0=log n-Log L, its analytic factors have coefficients

P0=1, P1=3L, P2=10L²-2Lb, P3=35L³-15L²b;

Q0=-1, Q1=(b-7L)/2, Q2=(-74L²+28Lb+b²)/6,
Q3=(-533L³+327L²b-9Lb²+b³)/12.

The box differentiates the WHOLE primitive, including log n. Its order-j coefficient is (j+1)Pj H0+(j+1)Qj-Pj. Each heavy-angle power b^k averages to t^k(k!)²/(2k+1)!.

Every light logarithmic moment is retained and reduced by the full integration-by-parts identity. With Mj=integral L^j, Jj=integral L^j Log L and a0=1-s/4,

Jj=[2j a0 J(j-1)-2Mj+2a0 M(j-1)]/(2j+1).

The boundary term vanishes because L=1 at the original endpoints. This is not a deletion of logarithmic moments. After reducing each distinct light channel and summing all ordered diagrams, every remaining J0 coefficient vanishes through n^-5.

Writing sigma2=s²+t²+u² and sigma3=stu, the complete base loop has

A1_base=g4/(16pi²)[F2/n²+F3/n³+F4/n4+F5/n5]+R6,

F2=-6 log n,
F3=-16 log n+10/3,
F4=(3 log n/2-157/90)sigma2-28 log n+10/9,
F5=(3 log n-13/21)sigma3+(12 log n-7969/630)sigma2-28 log n-2951/63.

Independent expanded and symmetric-polynomial derivations agree. All higher inverse-mass terms remain in R6.

## Direct joint holomorphy and complete product remainders

Use v=s+t/2-2, so s=2+v-t/2 and u=2-v-t/2. On the closed joint bidisk |v|,|t|<=1/4, the light parameter obeys |L-1|<=19/32 and Re L>=13/32>0; |L|,|b|<2. The full original parameter denominators have positive real parts. Thus the entire first-loop amplitude is jointly holomorphic on a neighborhood of this bidisk.

The S236 analytic-factor argument is explicitly extended to COMPLEX L,b with K=2. On the w disk of radius1/32, |P|<2 and |Q|<4. The nonanalytic -log w remains outside this Cauchy argument. Complete degree-three tails, weighted derivative tails and the exact rational vertex remainder give

|R_C|<=2^22 M/n5,
|R_D|<=3*2^23 M/n6,

with M=500 bounding the complete light logarithm factors at the actual hierarchy. The proof keeps every vertex/triangle/box cross-product, including the high products of retained polynomials. Per-channel constants are below2^27 M; their full sum is below10^12. Consequently

|R6|<R0=10^12 g4/(16pi²n6)

on the ENTIRE joint bidisk, not only sampled points.

## Three complete first-loop coefficient intervals

The coefficient conventions are b20=[v²t0], b21=[v²t1], b40=[v4t0], with derivative divisors2,2,24. Cauchy coefficient factors for R6 are16,64,256. The unchanged OS4 contact is constant and contributes zero to each observable.

The full first-loop shifts therefore satisfy

delta b20=g4/(16pi²)[(3 log n-157/45)/n4+(24 log n-7969/315)/n5]+e20,
delta b21=g4/(16pi²)(-3 log n+13/21)/n5+e21,
delta b40=e40,

where |e20|<16R0, |e21|<64R0 and |e40|<256R0. Exact elementary exponential estimates give394<log n<462 and retain the complete error margins.

At the actual parameters,

0<delta b20/(4lambda)<10^-203,
0<delta b21/(-3gamma)<10^-203,
|delta b40|/(gamma²/lambda)<10^-192.

The first b20 shift is positive and the first b21 shift negative. No sign is assigned to the first b40 shift. The tree-plus-first-loop b40 stays positive and within the stated relative interval around its separate tree value. These are explicitly TRUNCATED-amplitude coefficients, not exact physical LSZ/dispersion coefficients or all-loop Gram saturation.

## Independent validation and completed publication gate

The frozen package has17 source inputs and20 report fields,60 named identities,60 scalar entries,36 gates,nine controls and273 rejected unsupported inputs. It adds one matching record for93 records; the nine original primitive statuses are unchanged.

Independent checks cover complex analytic P/Q derivatives, full complex logarithmic moments, exact parameter moments, every product norm and both complex boundaries. Actual complete primitive tails are checked at1000 digits. An independent unexpanded full-loop numerical differentiation computes all three coefficients using85 digits and40-point parameter quadratures, with the correct coefficient factorials. All agree with the complete written bounds; the quadrature is diagnostic rather than an interval certificate.

The first private research probe used a structural comparison between factored and expanded forms; replacing it by exact rational cancellation passed. During initial implementation, a log1p-series helper was incorrectly given (1+d)/2 instead of the required increment (d-1)/2. The independent coefficient check immediately detected a Q0 mismatch of5/6. Correcting this private transcription and adding a dedicated helper-origin regression restored the intended formula BEFORE freezing. No frozen input changed.

Corrected private algebra passed, and initial science passed410 tests in1.71 seconds. The added full unexpanded-loop derivative test passed separately in31.70 seconds. Lint/format corrections, including explicit lambda-index binding, were completed before freezing. A private research note was moved outside the17-source manifest. Final preflight verified all sources, fields and counts in0.41 seconds and passed411 science tests in31.58 seconds. Repository science passed411 in31.83 seconds.

The46363-character native report was transferred in four checked chunks. All17 source hashes and every certificate field were independently checked. Report SHA:

`6d76badc63d988c9fae045c5a3628b20c7fb476ea9f45fdf98fd52a7a1b0e9fb`.

Fresh ordinary original-SymPy replay passed436 tests in2444.69 seconds. Standalone CLI replay passed. Full P8 regression passed49832 tests in4513.19 seconds, exit code0. The727-file snapshot SHA is

`b9bcf3673f8f24a8aa241369412c4b1212fd3c00f47cdace0d77417ca348d14c`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All128 original tuple self-checks passed. Final counters were36143 domain fallbacks,7340 exact descents and94 mixed fallbacks; collection used549 static namespace ancestors. Native, direct, ordinary and CLI retained original SymPy. Exact algebra and written complex-analysis proofs are not FORMALIZED.

Publication is limited to the exact21-file manifest, with every frozen source and report hash rechecked. Successors and unrelated P4/P9 work are excluded. New covariant-parent matching, omitted loops, full quantum UV and finite-gravity Regge control, the heavy state and same-state nonlinear bounce remain open. No user intervention is required.
