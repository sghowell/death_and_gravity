# P8 continuation: isolated shear pole, continuum and causal inverse

S6.223 completes the isolated homogeneous unit-Frobenius shear reference factor, including its pole and its ordinary finite-window causal inverse. It does not invert the full curved quantum system or establish a physical instability. Original P8 remains OPEN. The [full quantum-force graph](assessment-2026-09-12-p8-affine-quantum-forced-constraints.md) remains the coupled interface to solve.

## Original normalization, not a fitted finite term

Before the fixed factor64pi^2, the original Weyl coefficient is-1/30. The quadratic unit-Frobenius shear action therefore has fourth-derivative coefficient-1/60 and its Hessian has coefficient-1/30. Combining this unchanged finite term with the actual S199 Proca cut gives

F2(p)=-1/30-(p/30) integral_0^1 y^2(30-20y^2+3y^4)/[4m^2+p(1-y^2)]dy.

The physical cut is multiplied by64pi^2 in this unit-shear channel. Independently multiplying the actual trace cut by768pi^2 reproduces the old trace radial factor and its finite coefficient-4. The differing channel norms are retained. A cut alone does not fix either finite polynomial.

Write A2=-F2. Its strict half-plane imaginary-part sign and real-axis monotonicity show that it has exactly one simple first-sheet zero p=-r m^2 below threshold, and no other first-sheet or cut-bank zero. Exact rational moment enclosures give

2899/5000 < r < 5799/10000,
16 < R/m^2 < 17, where R=1/A2'(-r m^2)>0.

The residue of1/F2 is NEGATIVE R. The approximate values r=0.5798011453 and R/m^2=16.1801168580 are illustrations, not the root proof.

## Retain both the pole and the positive continuum

For tau>4m^2 let z=1-4m^2/tau and

D=-172/225+19z/30-z^2/10
  +sqrt(z)(30-20z+3z^2)atanh(sqrt(z))/30,
U=sqrt(z)(30-20z+3z^2)/60,
rho=U/(D^2+pi^2 U^2)>0.

The exact reciprocal representation is

1/F2(p)=-R/(p+r m^2)-integral_(4m^2)^infinity rho(tau)/(p+tau)dtau.

The positive static moments INCLUDE the pole:

R/(r m^2)+integral rho/tau=30,
R/(r^2 m^4)+integral rho/tau^2=675/(14m^2).

Threshold, ultraviolet and derivative estimates justify the unsubtracted contour and subsequent time-domain construction. In particular rho/sqrt(z) tends to50625/59168 at threshold, while rho log^2(tau/m^2) tends to60/13 at infinity. The explicit high-log tail bound rho<=50/u^2 for u=log(tau/(4m^2))>=16 is retained.

## Ordinary finite-window inverse, not a half-line L1 kernel

The inverse kernel is the negative sum of the undamped pole sinusoid and the continuum sine transform:

K_m(t)=-theta(t)[R sin(m sqrt(r)t)/(m sqrt(r))
                 +2 integral_(2m)^infinity rho(Omega^2)sin(Omega t)dOmega].

The continuum is an ordinary oscillatory integral. Its small-time bound is O(1/[t log^2(1/t)]) and its large-time bound is O(t^-3/2), after mass scaling. Thus the continuum is half-line L1. The complete kernel is L1 on every finite window, but NOT on the half-line because the nonzero pole is retained.

The absolutely convergent primitive J has-60<=J<=0 and J(0)=0. For prepared C1 input f with f(0)=0, K*f=J*f' and ||K*f||<=60T||f'||. General input retains the J(t)f(0) boundary term. No mass-independent half-line norm is inferred from scaling.

The forward causal distribution is defined without splitting divergent local/cut pieces:

F=-delta/30-D^2 G_m,
G_m=theta integral_0^1 W2(y)sin(Omega_y t)/[(1-y^2)Omega_y]dy,
Omega_y=2m/sqrt(1-y^2), ||G_m||<=9pi/(128m).

Laplace transformation of the absolutely convergent primitive proves both F*K=delta and K*F=delta. The finite-window graph includes the entire initial boundary; no initial atom or pole is removed. This is a complete inverse for the stated isolated factor, not a claim about poles of the full coupled operator.

## Independent checks and completed fresh validation

Independent finite radial moments enclose the root and residue without using a numerical root finder. A separate65-digit reconstruction checks four complex radial/closed-form values, four complete pole-plus-cut dispersions and both static moments. Omitting the pole produces a detected order-one error. Written threshold/tail and causal graph proofs accompany these diagnostics.

One initial private kernel check had an unfactored symbolic exact-zero expression; factoring it resolved the test without changing a formula. Final private science passed243 tests in4.09 seconds; repository science passed243 in4.14 seconds. Preflight verified17 inputs and20 fields in1.80 seconds. Formatting and lint passed.

Fresh original-SymPy ordinary replay passed268 tests in2316.27 seconds. Standalone CLI replay passed. Full P8 regression passed45203 tests in4475.92 seconds, final exit code0. Its699-file snapshot SHA is
`58a42c39eeea781300753ae3dd20aba52a355080f0be8085062ebca2c646fcf8`.
The full launcher retained the S219 helper-directory PYTHONPATH allowance. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36090 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The36173-character native report was transferred in four checked chunks. Its20 fields contain37 named identities,37 scalar entries,26 gates,nine controls and165 rejected inputs. Every frozen source/proof/test hash and the report SHA were checked:
`c7c0eb49287a4456328fe86b20a8362b70f540309439ba98b74ad1f191e0c1e6`.

The exact21-file publication manifest contains17 inputs, report, this audit, CLAIMS and README. The pending flat quotient checkpoint, private curved-reference work and unrelated P4/P9 changes are excluded. Functional and contour arguments are written proofs, not FORMALIZED.

## Remaining work

The next flat two-channel quotient checkpoint is still under fresh validation. The actual curved/state/contact/tree/matter normal form, compatible full quantum inverse, finite-coupling/nonlinear parent remainder, quantum background/stability, heavy/physical cutoff and original V/G/B remain open. All nine original primitive rows and prior matching rows are unchanged. Continuation remains active with no user-intervention blocker.
