# P8 continuation: original scalar unitarity and matching-scale tradeoff

S6.231 derives the complete original massive scalar tree in normalized identical channels and proves a conditional full-amplitude matching tradeoff. It does not establish a physical cutoff or exclude the original parent. Original P8 remains OPEN.

## Original target and identical normalization

The source-pinned S177 complete canonical functions and S182 vacuum are unchanged. The scalar mass is1, lambda=10^-600 and gamma=1024*10^-800. The Proca mass1000 and gravitational normalization sqrt(kappa)=10^400 are distinct scales. The nonconstant vacuum retuning jets below order2048 vanish; the constant vacuum term is retained separately.

For s+t+u=4, the COMPLETE tree is

A=2lambda[(s-2)²+(t-2)²+(u-2)²]+3gamma stu-8gamma.

The potential term is not discarded. The original vector-source and higher-field vertices do not contribute a connected four-scalar tree; this does not remove their real quantum corrections.

For the fully labelled amplitude, the full-sphere two-body phase space and the single identical-final-state factorial give

Im A(s,0)>=beta/(64pi) integral_-1^1 abs(A(s,x))² dx,
beta=sqrt(1-4/s).

Both identical incoming and outgoing state factors are retained. The normalized even partial waves are

t_l=beta/(64pi) integral A P_l,
A=(32pi/beta) sum_even(2l+1)t_l P_l,
S_l=1+2i t_l.

An exact unitary full channel matrix gives abs(S_l,elastic)<=1 even with inelastic channels, hence Im t_l>=abs(t_l)² and abs(Re t_l)<=1/2. An independent constant-vertex bubble and a coupled-channel matrix check fix the factor of two. A distinguishable16pi projection cannot be mixed with this labelled-identical optical coefficient. The [primary scalar-unitarity treatment](https://link.springer.com/article/10.1140/epjc/s10052-018-6127-z) supplies normalization context; the required inequality sign is derived directly here from the unit disk.

Writing A=a+b x² and Q0=a+b/3 gives t0=beta Q0/(32pi), t2=beta b/(240pi) and

rho_first=beta[Q0²+4b²/45]/(32pi).

This is the first elastic absorptive coefficient determined by the complete tree, not the full real one-loop amplitude.

## Nominal tree ceiling, not a physical cutoff

Exact shifted-polynomial inequalities show Q0>0, t0>=5abs(t2), and strictly increasing t0 from zero to infinity above threshold. Its unique real-part ceiling t0=1/2 lies at COM energy

10^133<E_U<2*10^133.

At energy10^134 the tree t0 exceeds50000. If an exact unitary amplitude exists there, its relative difference from this tree partial wave must exceed1-10^-5, by the exact unit-disk bound alone.

A nonzero real tree already lies outside the EXACT unitarity disk at lower energies: its imaginary part vanishes but its optical right side does not. E_U is therefore a nominal real-part ceiling, not the first nonunitary energy or a verified physical cutoff.

## Conditional full-amplitude dispersion implication

The physical premises are explicit and unproved for this original parent: an actual nongravitational S matrix with physical scalar pole mass1 and positive canonical LSZ normalization; the stated crossing and analytic neighborhood with correct known-light-pole accounting and no unresolved cut through it; and a twice-subtracted forward dispersion relation with a controlled vanishing arc. Unknown additional positive spectral weight is not erased.

In addition, suppose the FULL angular amplitude differs from the complete tree in L2 norm by at most eta times the tree norm at every s in[E²/2,E²], with eta<1. This error includes all light/heavy and higher quantum contributions. It is not deduced from the small tree coupling. Let physical b2=B''(0)/2 at the pole-subtracted crossing point and, optionally, require abs(b2-4lambda)<=4lambda delta.

The exact optical inequality, reverse triangle inequality and crossing kernel give

b2>=45(1-eta)² gamma² E^8/(131072pi²).

The central-angle lower bound used in the proof retains the full scalar mass and potential. At the actual parameters, pi²<10 implies

1+delta>9(1-eta)²(E/10^125)^8.

At E=10^125, eta<=1/2 and delta<=1 are incompatible with the physical premises: the lower bound gives b2>9lambda while the matching tolerance gives b2<=8lambda. Which premise or tolerance fails is not selected. This is a quantified necessary tradeoff, not a full UV-parent no-go.

The tree is uniformly smaller than10^-46 throughout the same physical angular range up to that energy. Thus this tradeoff is distinct from the nominal large-tree ceiling. Small tree size nevertheless does not prove small physical matching errors. The [improved-positivity analysis](https://arxiv.org/abs/1710.09611) provides primary context for this order-sensitive distinction; the massive normalization and actual constants are independently derived here.

## Relationship to the tensor and bounce frontiers

S230's retained tensor-pole frequency exceeds10^398 in the same anchored units. The hierarchy does not prove or exclude validity of that mean equation at its poles, delete a pole, choose new branch data or show a scalar repair of the tensor response.

Full quantum vacuum matching, finite-gravity Regge control and a common-parent nonlinear bounce remain open. A finite-band curved inverse or an all-momentum flat inverse is not a stability or unrestricted nonlinear theorem. No P8(a) scope is changed.

## Independent validation and completed publication gate

The package has17 frozen scientific inputs and20 report fields. Independent checks reconstruct the literal original massive amplitude, all phase-space and identical-state factors, exact angular projections, the first elastic cut, monotonicity, rational ceiling brackets and the full-error dispersion constants. It contains34 named identities,34 scalar entries,21 gates,nine controls and223 rejected unsupported inputs, with87 matching records and the same nine primitive statuses.

Before freezing, one private test compared algebraically identical expanded and factored crossing denominators structurally. Replacing that comparison by exact cancellation gave307 science tests passing in1.76 seconds; no amplitude or dispersion formula changed. The explicit physical mass/LSZ and crossing-cut premises, the nonconstant-retuning qualification and the primary2026 paper title were corrected before freezing. Lint/format passed. Preflight verified17 sources and20 fields in0.73 seconds and passed307 science tests in0.58 seconds. Repository science passed307 in1.99 seconds.

The native original-SymPy rebuild transferred the52192-character report in five checked chunks. Every source hash and the report hash are verified. Report SHA:

`34da90ca48a00689b3b9d7db12b367db703c970236caabcecd293b1061559cd5`.

Fresh ordinary replay with original SymPy passed332 tests in2483.98 seconds. Standalone CLI replay passed. Full P8 regression passed47527 tests in4549.12 seconds, exit code0. The715-file test snapshot SHA is

`73cda61005e9fb42db9a02b10c3bb884e112c36e1612d5b8a4012242f9af7645`.

Only full regression used the audited exact-GCD adapter and the frozen S219 helper-directory allowance. All128 original tuple self-checks passed; final counters were36151 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy. The written continuum proofs are not FORMALIZED.

Publication is restricted to the exact21-file manifest:17 frozen sources, report, this assessment, CLAIMS and root README. Frozen successors and unrelated P4/P9 work are excluded. No user intervention is currently required.
