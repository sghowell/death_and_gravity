# P8 continuation: scalar transfer and required spectral moments

S6.232 quantifies the full spectral weight and higher forward coefficient required by the original scalar matching tolerances. It retains all low cuts and unknown heavy positive weight. Its implications are conditional on explicit unproved physical analyticity and unitarity premises. Original P8 remains OPEN.

## Original coefficients and complete first elastic transfer

For the unchanged S177/S182 massive scalar target, use v=s+t/2-2 and the physical known-light-pole-subtracted amplitude B(v,t). Define b20=partial_v²B/2, b21=partial_t partial_v²B/2 and b40=partial_v4B/24 at v=t=0. The original TREE values are4lambda, -3gamma and0, where lambda=10^-600, gamma=1024*10^-800 and the scalar mass is1.

The original identical normalization from S231 is retained. The complete nonforward first elastic cut is obtained from a two-angle convolution, not the square of one nonforward amplitude. For A_tree=a+b x² it gives

rho_first(s,t)=(32pi/beta)[t0_tree²+5t2_tree² P2(1+2t/(s-4))],
partial_t rho_first(s,0)=beta b²/[60pi(s-4)].

The apparent denominator is regular at threshold because b contains(s-4)². Independent angular integration verifies this expression. At this first crossing-even jet, fixed-s and fixed-v derivatives agree; claiming otherwise would be a false negative control. Omitting the crossing half-shift in the dispersive denominator is a genuine error and changes the inverse-fourth moment.

## Explicit full-amplitude premises

The conditional argument requires an actual nongravitational S matrix with physical scalar mass1 and positive canonical LSZ normalization; no additional unresolved coupled s/u cut below4 after the specified known-light-pole accounting; fixed-t crossing and analyticity near t0; convergence of the positive partial-wave endpoint derivative; and a differentiated twice-subtracted dispersion relation with a controlled vanishing infinity arc.

These assumptions are stronger than a forward optical identity and are not proved for the original parent. Unknown heavy poles are retained as positive weight. A new low threshold or uncontrolled differentiated arc cannot be ignored. The [primary massive-Galileon analysis](https://arxiv.org/abs/1702.08577) supplies context for transfer positivity and additional operators; this checkpoint derives its own normalization and massive constants.

Writing rho=Im A_exact(s,0), rho_t=partial_t Im A_exact(s,0) and w=s-2, the conditional identities are

b20=(2/pi) integral rho/w³,
b21=(2/pi) integral[rho_t/w³-(3/2)rho/w4],
b40=(2/pi) integral rho/w5.

All cuts and retained atoms are included, and rho,rho_t>=0 under the hypotheses.

## Full low-weight requirement without erasing the light cut

Split at K=M²>4 into[4,K] and(K,infinity), without double-counting an atom at K. With J4_low=(2/pi) integral_low rho/w4, the complete positive decomposition gives

b21+3b20/[2(K-2)]>=-3J4_low/2.

If abs(b20-4lambda)<=4lambda delta0 and abs(b21+3gamma)<=3gamma delta1, with delta0>=0 and0<=delta1<1, then

J4_low>=2gamma(1-delta1)-4lambda(1+delta0)/(K-2).

At M=10^99, delta0<=1 and delta1<=1/2, the right side exceeds gamma/5. K is an analysis split, not a declared new threshold or physical cutoff.

An independent all-angle bound on the complete original first elastic contribution gives

0<J4_first<[48lambda²K+(9/16)gamma²K³+gamma²]/pi².

Every positive endpoint remainder in this upper bound is accounted for. At K=10^198 it is below gamma*10^-200. Thus the stated physical matching premises require a full moment more than2*10^199 times that original first elastic contribution. The original tree remains uniformly below10^-202 up to this split. No bound on the FULL low cut is assumed, and no claim that the required enhancement is impossible follows.

## A mandatory higher coefficient, with a positive-measure control

For the full positive measure dnu=(2/pi)rho ds/w³, the Gram matrix of1 and1/w gives J4_total²<=b20*b40. Since the negative b21 requires J4_total>=-2b21/3,

b40>=gamma²(1-delta1)²/[lambda(1+delta0)].

The named tolerances imply b40>=gamma²/(8lambda)>0. Exact matching of the two low coefficients would require at least gamma²/lambda. This is a required physical higher coefficient, not a computed loop coefficient.

A separately labelled rational positive atom shows compatibility of these moments at the level of positive measures. Set D=2lambda/gamma, M_H²=D+2 and g²=gamma D4=1/2^26. Its crossing-symmetric exchange function has b20=4lambda, b21=-3gamma, b40=gamma²/lambda and J4_total=2gamma. Its mass lies below the named split and it saturates the Gram inequality.

The atom does NOT reproduce the full original amplitude or independent functions. The real exchange tree lacks the continuous elastic absorptive part required by exact unitarity. It is not an exact S matrix, a width approximation, a common-parent bounce or an attained full quantum positivity optimum. The subsequent separate local model remains a distinct matching proposal.

## Independent validation and completed publication gate

The frozen package has17 scientific inputs and20 report fields,42 named identities,42 scalar entries,24 gates,nine controls and225 rejected unsupported inputs. It adds one matching record for88 records and leaves all nine original primitive statuses unchanged.

Private checks independently derive the nonforward cut, crossing kernel and derivative factors; test full low/high decompositions on positive multi-atom measures; verify the conservative massive integral and actual rational margins; and derive both the Gram bound and the distinct higher pole coefficient. The integral and conditional continuum arguments are written proofs, not FORMALIZED.

Before freezing, two private tests needed corrections to their checks: an expanded/factored denominator was compared by exact cancellation, and a scope assertion was aligned with the literal qualification. Corrected initial science passed310 tests in1.45 seconds. The added full Gram checks then passed323 in1.32 seconds. These changes did not relax the mathematical boundary. Lint/format passed. Preflight verified17 sources and20 fields in0.36 seconds and passed323 science tests in0.55 seconds. Repository science passed323 in1.66 seconds.

The63048-character native report was transferred in six checked chunks. Its SHA is

`6a6c17275c0ed6c2dec6e65deb596fc264945168d2ebd058a53277a13d0057ff`.

Fresh ordinary original-SymPy replay passed348 tests in2470.47 seconds. Standalone CLI replay passed. Full P8 regression passed47875 tests in4461.09 seconds, exit code0. The717-file snapshot SHA is

`a828c8942ca14b3e59951266e50211e0e079b14ef1dc14d2c3100d960b4954da`.

Only full regression used the audited exact-GCD adapter and frozen S219 helper-directory allowance. All128 original tuple self-checks passed; final counters were36133 domain fallbacks,7340 exact descents and94 mixed fallbacks. Native, direct, ordinary and CLI retained original SymPy.

Publication rechecks all17 source hashes and the report hash and is restricted to the exact21-file manifest. Frozen successors and unrelated P4/P9 work are excluded. Full quantum matching, finite-gravity Regge control, the nonlinear common-parent bounce and original V/G/B/P8 remain open. No user intervention is currently required.
