# P8 continuation: complete first-loop physical-angle remainder

S6.236 proves a uniform real-and-imaginary first-loop remainder for the unchanged separate V2S-T1-OS4 prescription. It retains every ordered diagram, the existing symmetric-value contact and the original massive target. It does not estimate omitted loops or replace the original affine/DHOST/Proca parent. Original P8 remains OPEN.

## A complete parameter primitive, not a momentum-integrand expansion

For L=1-s xi(1-xi) and b=t eta(1-eta), define

T(n,L,b)=integral_0^1 (1-z)/[nz+(1-z)²L-bz²] dz.

The complete triangle is the xi integral of T(n,L,0); the ordered box is minus the heavy-mass derivative of the xi,eta integral of T(n,L,b). Differentiating the entire expression retains the derivative of log n and both ordered channels.

The exact factorization Q=(L+alpha z)(1+beta z), with
delta=sqrt(n²-4nL+4Lb), alpha=(n-2L+delta)/2 and beta=(L-b)/alpha, gives a closed logarithmic primitive. Its light logarithm is Log(L-i0). The apparent beta=0 quotient has its analytic removable value; the branch sign is not changed.

Writing w=1/n gives T=w[P(w)(-log w-Log(L-i0))+Q(w)], where P and Q are analytic at zero. The nonanalytic -log w is retained explicitly, not subjected to a false Cauchy argument. For |L|,|b|<=K, on |w|<=1/(16K), the full analytic functions obey |P|<2 and |Q|<4.

The complete first two terms are

T=(H0-1)/n+[3L H0+(b-7L)/2]/n²+R_T,
-partial_n T=(H0-2)/n²+[6L H0-10L+b]/n³+R_D,

where H0=log n-Log(L-i0). If S>=4, K=S/4 and S/n<=1/8, complete Cauchy tails give

|R_T|<=64S²(|H0|+2)/n³,
|R_D|<=256S²(|H0|+3)/n4.

The light-threshold logarithm is integrable, including coincident roots. Factoring L=S(xi-r_plus)(xi-r_minus) bounds its full integral; no absolute-value integration through a Feynman double pole is used.

## Cancellation before taking bounds

With I(s)=1-s/6 and J(s)=integral L Log(L-i0), the triangle and box expansions retain all B and J terms. The heavy vertex factor is expanded with its exact rational remainder. Combining every bubble, triangle and ordered box cancels the complete light logarithmic terms through the first two inverse-mass orders.

The full three-channel result is

A1_base=g4/(16pi²)[-6 log n/n²+(-16 log n+10/3)/n³]+R_base,

with |R_base|<10^7 g4 S²/(16pi²n4). The displayed leading terms are independent of the on-shell kinematics; the bound includes every vertex cross-product and higher analytic tail.

The same decomposition at the symmetric point uses the auxiliary value S=4. Therefore the ALREADY FIXED OS4 contact cancels both universal terms, with no further finite condition. Throughout 4<=S<=10^196 and every real scattering angle,

|A1_OS4|<2*10^7 g4 S²/(16pi²n4),
|A1_OS4|/A_original<10^7 g²/(9n)<10^-199.

The original complete massive tree is strictly positive in this comparison and exceeds lambda(S-2)². Combined with S233, the tree-plus-first-loop expression differs from the original tree by less than 1/60+10^-199<1/59. This is explicitly a TRUNCATED-amplitude comparison, not an all-loop error or exact unitary S matrix.

A real physical-angle estimate alone does not bound analytic derivative coefficients at a different point. Those coefficients require a separate complex-neighborhood argument in a successor.

## Independent validation and completed publication gate

The frozen package has 17 scientific inputs and 20 report fields, 47 named identities, 47 scalar entries, 30 gates, nine controls and 243 rejected unsupported inputs. It adds one first-loop matching record for 92 records. All nine original primitive statuses remain unchanged.

Independent checks cover the exact primitive, full ordered box and heavy-mass derivative; positive and negative light parameters; the Feynman delta-function sign; complex-circle analytic bounds; complete threshold-root logarithms; and actual-hierarchy tails at 900 digits. A separate full unexpanded forward-loop quadrature reproduces the full tree-squared cut and checks its real remainder at two heavy masses. These diagnostics do not replace the continuum proofs.

Private algebra and cancellation probes passed. Initial private science passed 368 tests in 1.96 seconds. The added full-loop quadrature test passed independently with 368 deselected. There was no failed scientific run. Lint and format corrections were completed before freezing. Final preflight verified all 17 sources, 20 fields and counts in 0.43 seconds, and passed 369 science tests in 1.86 seconds. Repository science passed 369 in 1.91 seconds.

The 43732-character native report was transferred in four checked chunks; all 17 source hashes and every certificate field were verified. Report SHA:

`7a5e56979d2a755b7df1fc26f8e6d3d6a4e1bf760d4aa20911b139c250d28cb8`.

Fresh ordinary original-SymPy replay passed 394 tests in 2471.90 seconds. Standalone CLI replay passed. Full P8 regression passed 49396 tests in 4563.17 seconds, exit code 0. The 725-file snapshot SHA is

`f211298d2026782de31a76d82a1d303b0c96d8d378c582f7d5efe93d53d87af7`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All 128 original tuple self-checks passed. Final counters were 36155 domain fallbacks, 7340 exact descents and 94 mixed fallbacks; collection used 547 static namespace ancestors. Native, direct, ordinary and CLI retained original SymPy. Written complex-analysis bounds and exact algebra are not FORMALIZED proofs.

Publication is limited to the exact 21-file manifest, with all frozen source/report hashes rechecked. Successors and P4/P9 changes are excluded. Omitted loops, exact physical coefficient matching, a controlled heavy resonance, physical UV and finite-gravity Regge control, and the common-parent nonlinear bounce remain open. No user intervention is required.
