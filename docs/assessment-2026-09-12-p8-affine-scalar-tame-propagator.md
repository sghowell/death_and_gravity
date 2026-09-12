# P8 continuation: uniform two-chart classical scalar propagator

S6.221 replaces the compact-momentum classical comparison bound with an all-momentum polynomial bound for the unchanged fixed-QG1 scalar coefficient system. It does not invert the nonlocal quantum constraints. Original P8 remains OPEN. The [regular scalar Hamiltonian and full infrared pullback](assessment-2026-09-12-p8-affine-reduced-scalar-hamiltonian.md) is the source-pinned predecessor; the [phase erratum](assessment-2026-09-12-p8-retarded-phase-erratum.md) remains in force.

## Actual principal coefficients and complete charts

The actual polynomial background gives
Theta=H-t/(1+t^2)^4,
E=1-3/[2(1+t^2)^3],
ell=1/[10(1+t^2)^6], w=-ell E.
Set C=Theta E'-E Theta'+H E Theta-Theta^2 and F=C-w^2/2. Exact positive polynomial coefficients prove F>1/100 and Jc-F>1/1000 on the unchanged unit slab, with1/100<Jc<100. The principal squared speeds are1 and F/Jc. Both valid chart characteristic determinants are checked; this is not an argument obtained only by multiplying a singular formula byTheta^2.

The older core gamma complementary-chart algebra is explicitly source-pinned and reused, not presented as new. The new calculation retains the complete finite-transfer S220 Hamiltonian, full fixed retuning A,Tcorr,Jc, matter mode, weighted canonical boundary and q'=-2Hq.

Near the bounce, |t|<=1/4 gives |E|>1/4. The complementary coordinate b_c=-pv/(2q) is used there. It is NOT the physical S220 shift divergence b=pv/2. Away from the bounce, |t|>=1/8 gives |Theta|>=3/10 and the original curvature chart is valid. The two overlap, and switches at t=+-3/16 require at most two changes of chart.

Only |P|>=100 uses the two second-order charts; then q>=4096. With z=1/q, the central normalized denominator
Delta=(2E+3Tcorr z)^2-(2z+9Az^2)(2Jc+E^2ell^2)
lies strictly between1/8 and2. The complete kinetic and principal gradient matrices in both charts lie between1e-8 I and1e4 I. All48 entries of K,K',G,G',the skew velocity matrix and the remaining lower-order matrix have explicit polynomial enclosures below1e18 after their positive denominators are cleared.

These bounds use the actual coefficient jets through2 and the existing five controlled stress jets. They do not invent a numerical thirteen-jet stress estimate. Every q-dependent term and time derivative is retained. A separate exact control finds a central-chart degeneracy at bare-bounce q=152/25 while the original phase Hamiltonian is regular, which is why the low-momentum region is not put in that chart.

## Energy patching and the continuum comparison

The full second-order energy is
Energy=(ydot^T K ydot+q y^T G y)/2.
Its derivative retains K',G',q',H and the complete lower-order matrix; the skew velocity part does no work. Independent differentiation checks the identity even with a nonsymmetric lower-order matrix and nonzero forcing.

The complete coefficient bounds give a deliberately conservative energy-growth exponent below1e28. Forward and inverse phase/energy conversions, the two chart switches and the regular low-momentum first-order system yield

||U(t,s;P)||<=exp(1e29)(1+|P|^2)^6.

For every real r this supplies a continuous classical coefficient-sector propagator from H^(r+12) initial data and L1_t H^(r+12) additive mathematical phase forcing to C_t H^r. Plancherel, Duhamel and dominated continuity are proved in the notes. Twelve spatial derivatives are lost; the enormous constant is not a smallness, stability or same-space quantum contraction estimate. It cannot be combined with the small Gaussian coefficient by ignoring the derivative mismatch or the exp(1e29) factor.

## Independent checks and fresh immutable validation

Fourteen independently entered canonical-flow/chart fixtures check the complete second-order Euler equation, including nonzero rational test retuning. Six separate finite-q stationary-action saddle solves reconstruct the chart matrices and two independent original phase vectors. A deleted weighted-boundary control changes the equations. Further coefficient, characteristic, energy and phase-reconstruction checks supplement the continuous proof.

The first independent suite passed29 tests in74.04 seconds. The full48-entry/energy module probe passed in58.42 seconds. Before freezing, the central forward-velocity display was strengthened to retain the explicit matter-velocity contribution, without changing its500 bound. The two valid chart characteristic checks were added before freezing. Integrated science passed302 tests. Preflight checked18 inputs and20 fields in66.94 seconds, followed by302 science tests in6.63 seconds. Final private science passed302 in75.11 seconds; repository science passed302 in75.60 seconds. Formatting and lint passed.

Fresh original-SymPy ordinary replay passed327 tests in2339.52 seconds. Standalone CLI replay passed. Full P8 regression passed44649 tests in4443.19 seconds with final exit code0. Its695-file snapshot SHA is
`e682d475c4340953ac8658e8cd4ccbada6098b1e4f3bdcca79a7aacf0aa2e2eb`.
The launcher explicitly retained the S219 helper-directory PYTHONPATH allowance; no frozen test was edited. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, and final counters were36076 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The60116-character native report was transferred in six checked chunks. Its20 fields contain45 named identities,126 scalar entries,51 proof gates,nine controls and176 rejected inputs. All18 frozen source/proof/test hashes and the report SHA were verified:
`7560a2ba6fcb470f784662046e2561d641bf4f6f0034104a158d977093075255`.
The exact22-file publication manifest includes those18 inputs, report, this audit, CLAIMS and README. Pending S222/S223, private mixed-reference work and unrelated P4/P9 changes are excluded. Analytic energy patching and functional-space arguments are written proofs, not FORMALIZED.

## Remaining work

The quantum-force interface successor has passed ordinary and CLI replay and is still under full regression. An isolated shear pole-plus-cut inverse is also under fresh validation. Neither changes this checkpoint's classical scope. A complete nonzero-transfer curved coupled inverse, compatible source/response spaces, finite-coupling/nonlinear sourced-parent remainder, quantum background/stability, heavy/physical cutoff and original V/G/B remain open. Scoped P8(a), A.20-A.23, all nine original primitive statuses and all previous matching rows are unchanged. Continuation remains active; there is no user-intervention blocker.
