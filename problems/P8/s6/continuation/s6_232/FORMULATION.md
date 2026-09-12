# S6.232: conditional spectral and higher-coefficient requirements

## Original input and coefficient convention

The source-pinned original S177 fixed-canonical four-scalar tree belongs to the same S182 retained vacuum. Its nonconstant retuning does not change low field jets. The scalar mass squared is1, lambda=10^-600 and gamma=1024*10^-800. The Proca mass1000 is not the scalar mass.

Set v=s+t/2-2, u=2-t/2-v, and let B(v,t) be the physical amplitude after the specified known-light-pole subtraction. Define

b20=partial_v² B(0,0)/2,
b21=partial_t partial_v² B(0,0)/2,
b40=partial_v^4 B(0,0)/24.

The original TREE values are4lambda, -3gamma and0. They are not identified with physical quantum coefficients without an error premise.

## Explicit additional physical hypotheses

The conditional statements suppose an actual nongravitational scalar S matrix with physical pole mass1 and positive canonical LSZ normalization. After the specified original light-pole accounting there is no additional unresolved coupled s/u cut below4. The necessary fixed-t crossing and analyticity hold near t0; the positive partial-wave endpoint derivative converges; and the twice-subtracted dispersion relation can be differentiated in t with the relevant infinity arc vanishing.

These hypotheses strengthen the forward-only information. They are NOT proved for the original parent. A new lower threshold, nonconvergent endpoint derivative or uncontrolled transfer derivative of the arc requires new analysis. Neither tree mass nor the forward optical identity alone proves them.

Known original light poles are subtracted as specified. Unknown heavy poles are not erased: they remain positive spectral atoms or contribute their explicitly retained dispersive terms. All inelastic/light/heavy continuum weight remains. No finite-gravity massless pole is dropped to manufacture these nongravitational hypotheses.

## Exact conditional relations

Let rho=Im A_exact(S,0), rho_t=partial_t Im A_exact(S,0), and w=S-2. Under the hypotheses, both rho and rho_t are nonnegative and

b20=(2/pi) integral rho/w³,
b21=(2/pi) integral[rho_t/w³-(3/2)rho/w^4],
b40=(2/pi) integral rho/w^5.

All integrals range from4 to infinity and include retained atoms. For K=M²>4 split the measure into[4,K] and(K,infinity). Define J4_low with weight(2/pi)w^-4. The exact positive decomposition in notes/split.md gives

b21+3b20/[2(K-2)] >= -3J4_low/2.

For physical coefficient errors delta0>=0 and0<=delta1<1,

abs(b20-4lambda)<=4lambda delta0,
abs(b21+3gamma)<=3gamma delta1,

this requires
J4_low>=2gamma(1-delta1)-4lambda(1+delta0)/(K-2).

At M10^99, delta0<=1 and delta1<=1/2, its right side is strictly above gamma/5.

The positive full Gram matrix of1 and1/w yields J4_total²<=b20*b40. Since b21<0 under these error premises, J4_total>=-2b21/3. Therefore

b40>=gamma²(1-delta1)²/[lambda(1+delta0)].

The named tolerances give gamma²/(8lambda)>0. No physical coefficient is freely retuned in deriving these requirements.

## Independent comparison and boundary

The complete first elastic coefficient from the original tree has a conservative massive bound
0<J4_first<[48lambda²K+(9/16)gamma²K³+gamma²]/pi².
At K10^198 this is below gamma*10^-200. This is an upper bound on THAT contribution, not on the full amplitude. Together with the conditional full lower bound, it requires enhancement greater than2*10^199.

The separate rational positive-pole diagnostic matches only the two low coefficients and saturates the Gram bound as a positive measure. It is not an exact unitary S matrix, a physical width approximation, the original parent or a full-function/common-bounce matching construction. Passing these moments does not supply a UV completion.

The17 source-pinned inputs have exact algebra, rational inequalities, independent controls and written conditional continuum proofs. They are not FORMALIZED. Native/direct/ordinary/CLI retain original SymPy; only complete regression uses the audited exact-GCD adapter. No original primitive status is closed.
