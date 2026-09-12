# P8 continuation: a separate local heavy-scalar tree matching model

S6.233 constructs the separately named V2S-T1 classical two-scalar model and proves a complete massive, all-angle tree matching bound. It does not replace the frozen affine parent. Original P8 remains OPEN.

## Literal local action and the entire classical potential

Keep the original scalar mass 1 and targets lambda=10^-600 and gamma=1024*10^-800. Define D=2lambda/gamma, n=M_H²=D+2, g²=gamma D4=1/2^26, g=1/8192, and C=-g²(3/D-2/D²). The new action is

L=(partial phi)²/2-phi²/2+(partial H)²/2-nH²/2+(g/2)H phi²+(C/24)phi4.

Both kinetic terms are canonical. The full potential, not just its Hessian or a truncated valley, is

V=phi²/2+n[H-g phi²/(2n)]²/2+q phi4,
q=g²(D-1)/[6D²(D+2)]>0.

It is coercive with a unique global minimum at the origin. For the actual D>2, a second exact completion gives V>=(phi²+H²)/6 globally. Negative controls show why the contact cannot be omitted and why D>1 matters. Small g²/n and -C are classical parameter facts, not quantum error estimates.

## Complete tree matching and a physical angular window

The literal vertices give A_new=C+g² sum_channels 1/(n-channel). For s+t+u=4, set xi=channel-2. The complete finite geometric identity is

A_new=A_original+R,
A_original=2lambda sum_channels(channel-2)²+3gamma stu-8gamma,
R=gamma sum_channels xi4/(D-xi).

This uses the full massive original tree, not a massless approximation or two selected coefficients. On physical scattering kinematics, 4<=S=s<n and t,u<=0. If r=(S-2)/D, then R>0, A_original>0 and

0<R/A_original<=6r²/(1-r).

For center-of-mass energy E<=10^98, r<32/625 and the full pointwise relative error is less than 6144/370625<1/60. The complete original angular function is retained throughout the window.

At forward t=0, v=s-2, the two low coefficients b20=4lambda and b21=-3gamma match exactly. The separate model has b40=gamma²/lambda, the positive-measure lower value from S232. A nonzero constant shift 16gamma/(D+2) remains. Thus coefficient matching is not equality of the entire original amplitude or independent off-shell functions.

## All angular channels and the complete first elastic coefficient

Set k=(S-4)/2, a=n+k, d=C+g²/(n-S), and L=ln[(a+k)/(a-k)]. Then

A(x)=d+2g²a/(a²-k²x²),
I1=integral A dx=2d+2g²L/k,
I2=integral A² dx
  =2d²+(4dg²/k+2g4/(ak))L+4g4/(a²-k²).

The exact primitives and removable threshold limits are checked independently. The uniformly convergent even angular expansion has infinitely many nonzero even partial waves, as shown using Rodrigues projections. It cannot be replaced by the original tree's l=0,2 truncation.

With the original identical-particle normalization, t0=beta I1/(64pi) and rho_first=beta I2/(64pi). The full positive pointwise tree comparison gives

rho_first,original<=rho_first,new<(61/60)²rho_first,original

above threshold throughout the stated window, with both coefficients zero at threshold. The relative increase is less than 121/3600. This is the first elastic absorptive coefficient, not a bound on the real one-loop part, a full quantum remainder or an exact unitary S matrix.

The window lies below the heavy exchange pole and the heavy-pair threshold. Phi parity forbids a mixed phi-H two-body channel. The heavy scalar can decay to two light scalars; its real tree pole is not an exact stable quantum spectral atom.

## Scope and research boundary

This is a separately specified local classical action with a bounded potential, a complete tree matching theorem and a normalized first elastic comparison. It is not an all-loop construction, controlled resonance, exact S matrix, finite-gravity Regge result or common-parent bounce. The original affine/DHOST functions, quantum state and counterterms are not edited. The higher coefficient required by the conditional dispersion argument is supplied at tree level in the separate model, without claiming its full quantum value.

## Independent validation and completed publication gate

The frozen package has 17 scientific inputs and 20 report fields, 42 named identities, 42 scalar entries, 23 proof gates, nine controls and 213 rejected unsupported inputs. It adds one separate-model matching record, giving 89 records, while all nine original primitive statuses remain unchanged.

Independent tests reconstruct the full kinetic matrix, potential completions, literal vertex signs, massive angular identity, forward coefficients, ten Rodrigues projection cases, rational angular integrals and complete first-cut comparison. Three 80-digit quadrature fixtures and actual-parameter 1000-digit evaluations resolve the large contact/exchange cancellations and tiny remainder. Numerical diagnostics do not replace the written continuum bounds; the result is not FORMALIZED.

Before freezing, private angular and core checks needed exact cancellation for structurally different but equivalent expressions. The kinetic check was strengthened from a determinant-only assertion to the full matrix, and the explicit global coercivity gate was added. These were private validation improvements, not edits to frozen results. Private science passed 309 tests in 1.44 seconds. Lint/format passed. Final preflight checked the 17-source/20-field manifests and all counts in 0.44 seconds and passed 309 tests in 0.57 seconds. Repository science passed 309 in 0.61 seconds.

The 45977-character native report was transferred in four checked chunks. Its SHA is

`31716f2b701c0120c328dd8a2506f677eb30edeca956d531c57975e38d8d25cb`.

Fresh ordinary original-SymPy replay passed 334 tests in 2425.14 seconds. Standalone CLI replay passed. Full P8 regression passed 48209 tests in 4511.19 seconds, exit code 0. The 719-file snapshot SHA is

`b74de76b8e1db57aedfc416d89e8894ef48efc9880ae7ea80fa922c40952930d`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All 128 original tuple self-checks passed; final counters were 36148 domain fallbacks, 7340 exact descents and 94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

Publication rechecks all 17 source hashes and the report hash and is restricted to the exact 21-file manifest. Successor work and unrelated P4/P9 changes are excluded. Full quantum matching, common-parent bounce control, finite-gravity Regge control and original V/G/B/P8 remain open. No user intervention is required at this checkpoint.
