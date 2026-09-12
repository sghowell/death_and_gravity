# P8 continuation: full quantum-forced scalar interface and ordered graph equation

S6.222 derives the complete force-dependent scalar constraints and their ordered causal feedback equation. It does not establish existence of the full quantum inverse. Original P8 remains OPEN. The [uniform classical scalar propagator](assessment-2026-09-12-p8-affine-scalar-tame-propagator.md) and [full scalar Hamiltonian/infrared pullback](assessment-2026-09-12-p8-affine-reduced-scalar-hamiltonian.md) are the source-pinned inputs. The phase erratum remains in force.

## Keep both forces and both original chart contacts

Write g=(gn,gzeta,gb) for the additional scalar metric forces, normalized by the output density1/(kappa a^3). The physical scalar coordinates are(n,zeta,b), with zeta=v+delta n and b the physical shift divergence. This b is not the complementary S221 chart coordinate-pv/(2q).

The original scalar clock map has zeta(N,v)=v-log[1+2delta(N^-2-1)]/4. Its first lapse jet isdelta and its second is4delta^2-3delta. The nonzero Gaussian one-current therefore supplies the additional nn contact3a^3 pressure(4delta^2-3delta), besides the full ADM chart contact already in S220. Both are retained exactly once. Their omission cannot be justified by cancellation of the combined classical-plus-quantum background one-point equation.

The complete normalized force saddle, before eliminating velocities or auxiliary fields, has determinant-8Jc and gives

b=(pv-3gb)/2,
n=(Lc-gn-delta gzeta)/(2Jc).

The shift force cancels from the lapse only after using the complete momentum and shift relations. Simply inserting the classical b=pv/2 or n=Lc/(2Jc) would miss actual quantum-force dependence.

For Z=(v,sigma,pv,ps), the full interface is

Z'=KZ+Fg, s=CZ+Daux g,
F=-Jcan C^T,
Daux=-u u^T/(2Jc)-(3/2)e_b e_b^T, u=(1,delta,0).

Daux is symmetric of rank2, with eigenvalues0,-(1+delta^2)/(2Jc),-3/2 and null direction(-delta,1,0). The interface reconstructs all four weighted phase equations. It is not a four-ODE quantum Hamiltonian: g retains the full memory and higher-derivative response.

## Exact ordered graph equivalence, not inverse existence

Let Rhat denote the complete scalar reference response with both contacts, and set Qbar=(kappa a^3)^-1 Rhat. This is an OUTPUT density normalization; multiplication by that time-dependent factor does not commute with the retarded response.

Using optional mathematical drives h in the phase equation and e in the scalar response, let G0 be the S221 zero-past classical propagator and Zh=G0h. Define

Taux=Daux+C G0F.

Classical variation of constants and the exact force interface give

[I-Qbar Taux]g=e+Qbar C Zh,
Z=Zh+G0Fg, s=CZh+Taux g.

An ordered block factorization proves both directions of this equivalence on the common graph domain. It does not assume Qbar, Daux or the Schur operator is invertible.

Precisely, use the initial-prepared source completion X_r=H13_t H^(r+8)_x and Y_r=L2_t H^(r-3)_x. The graph consists of g inY_r for which Taux g belongs toX_r, with e inY_r andC Zh inX_r. Smooth compact-momentum phase drives with the original initial zero germ give one explicit admissible preparation class. No density, invariance, surjectivity or bounded inverse of this graph is claimed.

The source-side current bound has no detector time derivative, so Riesz representation and spatial translation give Qbar:X_r toY_r with bound1e-683. This is a completed weak-response mapping, not a source-space self-map.

## What the current first-response bounds do and do not imply

The complete force map satisfies ||F||<400(1+P^2), and ||Daux||<100. The classical propagator therefore gives

G0F:L2H^s toC H^(s-14), bound400 exp(1e29),
Taux:L2H^s toL2H^(s-16), bound2e5 exp(1e29).

For a prescribed prepared scalar source S0 inX_r, its first force g1=Qbar S0 lies inY_r. Its classical phase response Z1=G0Fg1 lies inC H^(r-17), bounded by exp(1e29)1e-680||S0||. The returned scalar source Taux g1 is controlled only inL2H^(r-19), with bound exp(1e29)1e-677||S0||. That is not the required H13H^(r+8) domain.

When S0 is itself the S220 comparison-phase source, the unchanged finite high-jet coefficient C13 is retained and the phase-to-phase estimate loses27 spatial derivatives. No numerical thirteen-jet stress bound is inferred from five known stress jets. The giant classical constant and derivative mismatch prevent a claimed small Neumann contraction from these estimates.

This gap is not an instability or nonexistence proof. An exact comparison multiplier has the same kind of derivative-losing mapping and fails the proposed same-space bound, while a different complete operator built from it still has a bounded inverse. Formal second-response coefficients are bookkeeping only; no finite-coupling Born remainder is asserted.

## Independent checks and completed fresh validation

Nine independent literal action saddles cover negative, zero and positive Theta and transfers1e-12,2,10000. They check both momenta, both forced constraints, all four weighted equations and the physical metric tuple. Three separate nonlinear clock Hessians detect omitted or duplicated one-current contacts. A two-time rational causal matrix fixture compares the full block solve with ordered Schur reconstruction and detects a deleted auxiliary block or reversed density order.

An initial private block fixture used Python floating division in one nominally exact coefficient. Converting it to an exact rational fixed the fixture without changing the physical formulas. Independent science passed18 tests in2.05 seconds. After adding the actual fixed phase-map bridge, the integrated science suite passed261 tests. Preflight verified17 inputs and20 fields in6.53 seconds; final private science passed261 in7.81 seconds and repository science261 in7.92 seconds. Formatting and lint passed.

Fresh original-SymPy ordinary replay passed286 tests in2307.78 seconds. Standalone CLI replay passed. Full P8 regression passed44935 tests in4478.66 seconds, final exit code0. Its697-file snapshot SHA is
`61da73e3f520e653f65edba3296fca6bea073edebc2060550cf368ce2c4e913b`.
The full launcher retained the explicit S219 helper-directory PYTHONPATH allowance. Only full regression used the audited exact-GCD adapter:128 original tuple comparisons passed, with final counters36085 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

The43512-character native report was transferred in four checked chunks. Its20 fields contain41 named identities,195 scalar entries,37 gates,nine controls and164 rejected inputs. Every frozen source/proof/test hash and report SHA was checked:
`05e36391c6636fa167c5569d79f706f04029548d9517c0686fa29aba70b15a40`.
The exact21-file publication manifest contains17 inputs, report, this audit, CLAIMS and README. Pending S223/S224, private curved-reference research and unrelated P4/P9 changes are excluded. The functional analysis is written proof, not FORMALIZED.

## Remaining work

The isolated shear pole-plus-cut inverse and full flat scalar quotient reference inverse successors are under fresh validation. They do not retroactively turn the present graph identity into a full curved coupled inverse. The complete curved/state/contact/tree/matter normal form, compatible-space inverse, finite-coupling/nonlinear sourced-parent remainder, quantum background/stability, heavy/physical cutoff and original V/G/B remain open. Scoped P8(a), A.20-A.23, all nine original primitive rows and prior matching rows are unchanged. Continuation remains active with no user-intervention blocker.
