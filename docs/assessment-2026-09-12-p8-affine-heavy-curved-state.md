# P8 continuation: specified heavy curved state, complete stress and fixed mean profile

S6.240 supplies an actual added-heavy Gaussian state on the curved clock reference and a full covariant scalar renormalization prescription. It proves all-momentum bounds for the complete reference stress and its first five time derivatives, then fixes a new mean-matching profile. This closes the specified added-heavy one-point input, not the complete quantum response, interacting clock, nonlinear bounce, UV, Regge or original P8.

## Actual state and preparation

For the existing heavy mass squared n=10^200/512+2 and physical a=(1+t²)², choose the state of low energy with the explicit positive normalized smooth sampling supported on[-3/4,-1/2]. Its modes solve the exact massive KG equation and have the physical annihilation convention S conjugate(S')-conjugate(S) S'=i/a³. The full Bogoliubov formula is phase-regular at vanishing off-diagonal sampling coefficient, has the correct CCR, minimizes the sampled energy and is independent of the initial normalized basis up to the irrelevant global phase.

The state is prepared once. Later compact metric/source variations retain the same covariance in the shared initial neighborhood, propagate it with the changed equation and retain the causal coherent mean; no live reminimization is allowed. [Olbermann's state-of-low-energy construction](https://arxiv.org/abs/0704.2986v2), Theorems3.1 and4.9, supplies the regularity/Hadamard input. The quantitative bounds are derived here, not attributed to that theorem. The sign convention is checked independently.

## Entire covariant scalar prescription

The complete scalar adiabatic orders0,2,4 and their energy/pressure Ward identities are derived in spatial dimension D=3-2epsilon. The radial integrals retain all D-dependent tensor and coefficient factors. The complete covariant pole action is varied BEFORE taking D3, including the finite evanescent contribution from that variation. Setting D3 too early or subtracting component poles separately changes the finite answer and is rejected.

All six finite energy/pressure components match the full scalar action
[(3/2-log n)n²+(log n-1)n R_old/3-2log n a2]/(64pi²),
where a2=R_old²/72+(Riemann²-Ricci²)/180. The different old vector finite coefficients are not borrowed. [The scalar adiabatic/DeWitt-Schwinger comparison](https://arxiv.org/abs/1412.7570v2) provides context; the exact dimension-dependent variation and finite components used here are independently computed.

The renormalized reference stress is the FULL exact-state mode integral minus its full order-four subtraction, plus that fixed finite local action. This is not a finite WKB stress substituted for the exact state.

## Complete uniform remainder and full momentum integral

A common complex time neighborhood gives a holomorphic frequency branch with explicit mass-and-momentum lower bound Omega_star²=n+p²/16. Seven full Riccati iterates, with explicit losses of time-analytic width, bound the complete sixth-iterate defect by10^60 Omega_star^-12. Exact variation of constants bounds the actual normalized mode error. Ten integrations by parts use the entire compact sampling function and every endpoint; the full SLE mixing satisfies |beta|<10^83 Omega_star^-11.

The energy-scaled exact oscillator system and both physical stress matrices give explicit derivative recursions through five. The full auxiliary readout is phase-cancelled before Cauchy estimates; an oscillatory complex phase is not bounded as if it were harmless. The entire order-four subtraction tail is controlled with the retained1/2048 time margin. Exact state, exact evolution and full auxiliary-tail errors combine to a10^301 Omega_star^-5 modewise bound. The complete radial integral is finite at every momentum and gives the uniform reference state/subtraction stress-jet bound10^310/n.

The complete finite local action is also bounded, including its n² vacuum term and all curvature derivatives. Together the normalized added-heavy reference stress and its first five derivatives are below10^-400 on[-1,1]. The source-defined integrals, absolute convergence and derivative interchange are retained.

## New fixed QG2 mean profile, not automatic stability transfer

The explicitly named CD-REG-AFFINE-ISO-QG2-H8A420 candidate uses a smooth clock profile fixed from those full reference integrals. It cancels the added-heavy mean stress and satisfies the clock Ward equation. The original Proca/M1 choices and S238 recomputed physical affine mean are retained.

The two free scalar Gaussian vacuum constants are cancelled in the vacuum region. Nonconstant profile terms begin at degree2048 there, so the complete S239 first-loop vacuum matching is unchanged. This does not prove the interacting quantum effective potential or the interchange of quantization with the gravitational limit.

A full complex-tube/Cauchy estimate, including all derivatives and rational localizer denominators, gives a NEW combined four-jet budget10^-399. The older classical10^-2700 allowance is not reused. R and the affine maps are unchanged, but the bare clock F Hessian changes. Its classical comparison hypotheses and the full quantum response require a fresh audit; no old stability or inverse result is transferred automatically.

## Independent validation and private corrections

Independent nonlinear lapse/scale Euler calculations check1,R,R²,Ricci²,Riemann² in spatial dimensions3,4,5. Independent arbitrary-D pointwise adiabatic Ward tests vanish at all three retained orders. High-precision exact SLE numerical fixtures check normalized bases, the phase-regular limit, quadrature convergence and Cauchy evolution. They are diagnostics; the actual enormous-mass and all-momentum conclusions use the written uniform inequalities.

Actual full seven-iterate WKB probes at1800 and2000 digits agree and retain a nonzero scaled sixth defect. A private auxiliary readout initially used inv_a^(3/2) where the physical density requires inv_a³; independent review corrected it before freezing and added direct physical energy/pressure tests plus a wrong-volume-factor negative control. The defect estimate was unaffected. The pressure auxiliary bound was strengthened to use81+192 rather than an energy-only81+64 estimate. The time-Cauchy radius was corrected to1/2048, and the inadequate1/4096 margin is retained as a negative control. Lint/import/unused bindings and formatting were resolved before freezing. No frozen source, proof, test or report was changed.

Final private preflight checked17 sources and20 report fields in0.747550 seconds; final private science passed527 tests in3.76 seconds. Repository audit completed in0.747509 seconds and science passed527 tests in3.85 seconds. The report has57 named scalar identities,74 proof gates,nine controls and343 rejected unsupported inputs. The nine original primitive rows and95 previous matching records remain unchanged; one separate state/mean record is appended.

The78208-character original-SymPy native report was transferred in seven checked chunks. Its exact17-source hashes, private/repository bytes and20-field manifest were verified. Report SHA:

`d4a6081d0482b2168255299cdd0080c90a1716eeffbc70f2645e6e2452e78eda`.

Fresh ordinary replay passed552 tests in2512.85 seconds and standalone CLI replay passed. The complete P8 regression passed51375 tests in4524.19 seconds, exit code0.

The full regression's733-test-file snapshot SHA is
`72ee233903221230e0ca383f9f6fd4f505dfc3f207a275f5b7c1a97f7cbd9dd5`.
All128 exact-GCD original-tuple self-checks passed; collection used555 static namespace ancestors. Final adapter counters were36154 domain fallbacks,7340 exact descents and94 mixed fallbacks.
Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. Native, direct, ordinary and CLI retained original SymPy. The algebra and written proofs are not FORMALIZED.

Publication is restricted to the exact21-file manifest, with every frozen source/report hash rechecked. The private S241 successor and unrelated P4/P9 changes are excluded. Full heavy quantum response and constraints, interacting light/mixed clock loops, quantum gravitational decoupling, nonlinear same-state bounce, physical UV, finite-gravity Regge and original V/G/B/P8 remain open. No user intervention is required.
