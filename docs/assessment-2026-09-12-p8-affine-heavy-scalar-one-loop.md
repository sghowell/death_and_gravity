# P8 continuation: exact heavy integration and controlled first-loop input

S6.234 advances the separately named V2S-T1 model from its classical tree to a complete finite-regulator heavy integration, explicit first light self-energy, complex-disk remainder and first heavy decay coefficient. It is not an original-parent replacement or full quantum matching. Original P8 remains OPEN.

## Exact finite-regulator action and complete light Hessian

For K_H=-Delta+M_H² and G_H=K_H^-1, the Gaussian heavy integration gives

S_eff=quadratic_phi-C integral phi4/24-g²<phi²,G_H phi²>/8.

The heavy determinant is independent of phi in flat space, but is metric dependent on a curved background. The operator bound G_H<=1/M_H² retains the entire positive classical quartic margin q from S233. This is an exact finite-regulator statement at fixed bare parameters, not a continuum quantum measure or a counterterm-uniform limit.

The complete light Hessian correction is

W=-C M_phi²/2-g² M_(G_H phi²)/2-g² M_phi G_H M_phi.

The original ordering remains. An independent full two-field Schur complement verifies it. The degree-two and degree-four first-loop field jets are Tr(G_phi W)/2 and -Tr(G_phi W G_phi W)/4. Mixed heavy/light loops remain despite the field-independent heavy determinant. The [primary functional matching treatment](https://arxiv.org/abs/1604.01019) provides context; the actual model and Hessian are derived here.

## Heavy one-point condition and explicit renormalization convention

The first Euclidean heavy one-point term is -g I_phi H/2. The counterterm j1 H with j1=g I_phi/2 cancels it. Eliminating H with source g phi²/2-j1 supplies the corresponding quadratic cross term, cancelling the local heavy-source light tadpole. The quartic tadpole and mixed bubble remain; deleting either is not the stated calculation.

The separate scheme uses dimensional regularization in4-2epsilon, MSbar UV subtraction at mu=1, zero heavy one-point function and finite light on-shell mass and residue conditions at s=1. Tree parameters remain those of S233. There is no finite four-point condition in this checkpoint.

Define the 1PI insertion as iPi, so the inverse convention is s-1+Pi. With the standard B0 notation,

Pi_MS(s)=-C/(32pi²)+g² B0_MS(s;1,M_H²)/(16pi²),
B0_MS=-integral_0^1 Log[x M_H²+(1-x)-x(1-x)s-i0]dx.

The [primary scalar-integral conventions](https://arxiv.org/abs/0709.1075) fix notation, while independent Wick routing fixes all factors and signs here. The mixed bubble has no identical-particle1/2; the heavy light-pair bubble does.

## Entire subtracted logarithm and a complex-disk estimate

Let F=(1-x)²+x M_H² and alpha=x(1-x)/F. The fully subtracted first-loop insertion is

Pi_OS(s)=g²/(16pi²) integral[-Log(1-alpha(s-1))-alpha(s-1)]dx.

Its zero value and slope at s=1 are imposed FIRST-ORDER renormalization conditions, not exact LSZ conclusions. Complete parameter bounds give |Pi_MS(1)|<10^-7 and0<Pi_MS'(1)<10^-207.

For |d=s-1|<=R and r=R/M_H²<1, alpha<=(1-x)/M_H² and the full logarithm series yield

|Pi_OS(s)|<=g²|d|²/[96pi² M_H4(1-r)].

At R=10^196, r<32/625. Writing the analytically extended q=Pi_OS/d, its uniform norm is below an explicit epsilon approximately4.7651*10^-210. The algebraic reciprocal of the ONE-LOOP-TRUNCATED inverse d(1+q) has only its anchor simple pole with unit residue in this disk, and its ratio to1/d differs by less than10^-209.

This is not the exact quantum propagator, an all-loop resummation or a four-point angular error estimate. The mixed bubble retains the cut[(M_H+1)²,infinity). Exact threshold-square identities show why the pseudothreshold is not a first-sheet cut of this bubble. They do not exclude lower higher-loop light cuts.

## Complete first absorptive coefficients and heavy instability

Independent discontinuity calculations give

Im Pi_H(s)=g² sqrt(1-4/s)/(32pi), s>4,
Im Pi_phi(s)=g² sqrt[(s-(M_H+1)²)(s-(M_H-1)²)]/(16pi s)

above the mixed threshold. Combining the first heavy coefficient with the correctly normalized phase space, flux and identical factor gives

Gamma_H,first=g² sqrt(1-4/M_H²)/(32pi M_H),
10^-208<Gamma_H,first/M_H<10^-207.

The formal outgoing second-sheet pole shift has negative imaginary part. No exact width, controlled full resonance or stable heavy quantum atom is claimed.

## Validation and completed publication gate

The frozen package has17 source inputs and20 report fields,58 named identities,58 scalar entries,25 gates,nine controls and228 rejected unsupported inputs. It adds one first-loop matching record for90 records; all nine original primitive statuses remain unchanged.

Independent finite-matrix and Gaussian calculations verify the full action, Hessian, Schur ordering and trace-log jets. Labelled Wick routes, source-counterterm cancellation, literal insertion signs, complex B0 quadratures, threshold squares and both discontinuities are checked separately. Actual-parameter diagnostics use500 decimal digits and resolve the heavy endpoint layer. The continuum bounds are written logarithm and quadratic-form arguments, not numerical samples or FORMALIZED proofs.

Before freezing, a private test file had one stray patch marker causing a collection error. Removing that marker restored the intended code; science passed354 tests in47.70 seconds. The logarithm bound was strengthened with its exact positive exponential-series check. Lint/format passed. Final preflight verified17 sources and20 fields in0.35 seconds and passed354 science tests in46.89 seconds. Repository science passed354 in47.14 seconds.

The47711-character native report was transferred in four checked chunks. Its SHA is

`d6646d25d158bcfeb32bf07895e8fcc2ca08e05988919f803ce50aa5e45f58a5`.

Fresh ordinary original-SymPy replay passed379 tests in2489.85 seconds. Standalone CLI replay passed. Full P8 regression passed48588 tests in4510.72 seconds, exit code0. The721-file snapshot SHA is

`e3075651efe1da0daa90899f891f87a744ed024cf956b5e129b9267f0aa6e4d8`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All128 original tuple self-checks passed; final counters were36138 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

Publication is restricted to the exact21-file manifest, with every frozen source and report hash rechecked. Successors and unrelated P4/P9 changes are excluded. Full quantum matching, finite-gravity Regge control, the nonlinear common-parent bounce and original V/G/B/P8 remain open. No user intervention is required.
