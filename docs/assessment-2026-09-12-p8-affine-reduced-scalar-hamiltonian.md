# P8 continuation: regular scalar comparison and full infrared pullback

S6.220 supplies the complete current fixed-QG1 scalar coefficient-sector Hamiltonian and the infrared-safe pullback of the full reference Gaussian metric response. Original P8 remains OPEN. The [complete ordered spatial response](assessment-2026-09-12-p8-affine-full-spatial-remainder.md) is its source-pinned predecessor. The [phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) remains in force; no withdrawn physical formula is restored.

## Full nonzero-transfer scalar action and current retuning

The literal ADM reduction retains the independent matter wave, lapse, shift divergence b, spatial curvature, scalar mixing and weighted time/spatial boundaries. It is not the homogeneous fixed-charge Routh problem. With q=|P|²/a², the full coefficient action divided by kappa a³ is
L2=-3vdot²+(J+w²/2-3Theta²)n²+6Theta n vdot
   +sigmadot²/2+w n sigmadot-3ell vdot sigma
   +2b(vdot-Theta n)+ell b sigma
   +qv²+2Eqnv-q sigma²/2.

The fixed S182 profile is differentiated as its complete order1024 switch, not replaced globally by a polynomial. Its local quadratic correction is
DeltaJ n²+3Tcorr nv+9Av²/2,
where A=-pressure_ref/kappa,
DeltaJ=(21delta²-3delta)A/2+(1-6delta)B,
B=-(rho_ref+pressure_ref)/(2kappa),
and Tcorr=(1+3delta)A-2B.
Keeping only the pivot shift would omit actual mixed and volume terms. The full profile remains fixed under metric variation.

The original positive J and the fixed five-jet stress bound give Jc=J+DeltaJ>1/100 on the unit slab. Legendre transformation before eliminating the two auxiliary variables yields
b=pv/2, n=Lc/(2Jc),
Lc=Theta pv-w ps+3Theta ell sigma-(2Eq+3Tcorr)v.
The complete reduced coefficient Hamiltonian is
Hc=ps²/2-ell pv sigma/2+(q/2-3ell²/4)sigma²-qv²-9Av²/2+Lc²/(4Jc).
It is regular at Theta=0. The normalized momenta obey p'+3Hp=-H_field. This coefficient-sector elimination is not elimination of the nonlocal quantum lapse/shift equations.

Independent saddle solves check positive, negative and zero Theta at very small and large q. Separate three-dimensional extrinsic-curvature spatial means and complete fixed-profile derivatives check the underlying action. At the bounce, q=1,v=ps=0,sigma=1,pv=20 gives Hc=-203/400, a control against claiming the original instantaneous Hamiltonian is positive. It is not a time-dependent instability proof.

## Infrared cancellation requires both ordered Ward identities and contacts

For scalar metric inputs(n,zeta,b), zeta=v+delta n, the shift itself contains inverse momentum. The full source-retarded and detector-advanced synchronous maps instead give
eta=I n, c=I b-|P|² I(a^-2 eta),
Qs=2(zeta-H eta)I-2c Pi,
with Pi=Phat Phat^T and operator/Frobenius norm1.

The complete ordered source-density and ADM-chart terms cancel the inverse-square shift pole. Source covariance and detector conservation are distinct identities; neither is replaced by symmetry of the retarded kernel. Removing the ADM chart leaves the nonzero term2a²A_density bD cG'/|P|². An independent coordinate-density push detects this divergence at momentum1e-15, while the complete form remains finite. A separate detector pull checks the other ordering.

The actual clock metric is Qphysical=2vI-log R(N)I/2. Its nonzero second lapse jet requires the additional one-current contact
3a³ pressure_ref (8delta²-6delta)nD nG/2.
This is additional to the full ADM chart, and cannot be deleted using the combined one-point cancellation.

The completed S219 response consequently obeys
|Rscalar(D,G)|<1e117 V03[D] U138[G].
The detector has three spatial derivatives and no time derivatives; the source has eight spatial derivatives and time jets through13. The proof first works with a Fourier hole, then extends uniquely by the hole-independent estimate and density. No unweighted shift norm or direction-independent projector limit is assumed. A measure-zero P=0 assignment does not solve the literal homogeneous constraint.

The spatial anchored approximation, together with exact Ward and clock terms, has a1e76/K scalar error. This is not a Ward identity for a finite sharp lapse/shift band, which is not diffeomorphism invariant.

## Phase comparison and remaining derivative mismatch

The actual comparison map is s=L0(t)Z+|P|²L2(t)Z. It obeys V03[sD]<=400 V05[ZD] and U138[sG]<=C13 U13,10[ZG], where C13 is the explicit finite sum of coefficient jets through13. Smoothness and Jc>1/100 establish finiteness; no numerical thirteen-jet stress bound is inferred from five controlled jets. Thus
|Rphase(D,G)|<1e120 C13 V05[D] U13,10[G].

The complete classical coefficient generator is bounded by1000(1+|P|²)^2. On a mathematical compact transfer band |P|<=Lambda, its propagator is bounded by exp[1000(1+Lambda²)^2|t-s|]. This gives a regular classical comparison through the bounce, but not a uniform continuum estimate, physical cutoff or quantum inverse. The derivative-losing weak response cannot be treated as a same-space contraction using this display.

## Independent checks and immutable fresh validation

Exploration corrected an expected saddle-pivot factor, a shell-launch quoting error, an unfactored exact-zero assertion and an expected arithmetic margin before freezing. No physical coefficient was changed to hide a failing control. Independent science first passed32 tests, then41 after expanded ADM/detector checks. Integrated science passed317 tests. Exact preflight checked17 inputs and20 fields, with317 science tests passing in5.38 seconds. Final private science passed317 in21.46 seconds; repository science passed317 in21.27 seconds. Formatting and lint passed.

Fresh original-SymPy ordinary replay passed342 tests in2325.24 seconds. CLI replay passed. Full P8 regression passed44322 tests in4333.13 seconds with final exit code0. Its693-file snapshot SHA is
`66846eb9a724fe1bfe53d103b0b7a48b05b8d9f1d8dfffede1cccbb189e12a77`.
The full launcher explicitly retained the S219 helper-directory PYTHONPATH allowance; frozen tests were not edited. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36057 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The52568-character native report was transferred in five checked chunks. Its20 fields contain62 named identities,115 scalar entries,55 gates,nine controls and157 rejected inputs. All17 frozen source/proof/test hashes and report SHA were verified:
`cfd123d4faed1ce2efddb229931f6af79df203a7ea3c5af091d9d0e1857db4d8`.
The exact21-file publication manifest includes those17 inputs, report, this audit, CLAIMS and README. Pending S221, private quantum-force work and unrelated P4/P9 changes are excluded. Continuous Sobolev, preparation, smooth-coefficient and propagator arguments are written proofs, not FORMALIZED.

## Remaining work

A uniform classical propagator successor is under fresh validation. The exact quantum-forced constraints and compatible-space coupled inverse remain separate tasks, as do the nonlinear sourced-parent remainder, quantum background/stability, heavy/physical cutoff and original V/G/B. Scoped P8(a), A.20–A.23, all nine primitive statuses and all previous matching rows remain unchanged. There is no user-intervention blocker; continuation is active.
