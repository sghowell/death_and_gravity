# Conservative massive bound for the ORIGINAL first elastic moment

This note bounds only the complete first elastic absorptive coefficient fixed by the original tree. It is NOT an upper bound on the full quantum cut. Additional light/heavy channels, higher loops, real matching effects and nonperturbative contributions have not been bounded.

For S>=4 and abs(x)<=1, the complete original tree obeys

abs(A_tree(S,x))<=B(S)=4lambda S²+(3gamma/4)S³+8gamma.

Indeed the lambda angular factor is positive and bounded above by4lambda S². The Galileon factor has abs(1-x²)<=1 and(S-4)²<=S², while the potential contributes its full8gamma in absolute value. No threshold or constant term is discarded in forming this majorant.

The exact identical-state optical factor from S231 gives

rho_first=beta/(64pi) integral_-1^1 A_tree² dx<=B(S)²/(32pi).

For the inverse-fourth moment, S-2>=S/2, so

J4_first(K)=(2/pi) integral_4^K rho_first/(S-2)^4 dS
          <=(1/pi²) integral_4^K B(S)²/S^4 dS.

Use(a+b+c)²<=3(a²+b²+c²). The difference is the sum of the three pairwise squared differences. Exact integration and positive discarded lower-end terms give

J4_first(K)<[48lambda²K+(9/16)gamma²K³+gamma²]/pi².

For example integral_4^infinity S^-4 dS=1/192 controls the full potential term. Before replacing the finite integrals by this upper expression, the positive difference is
[192lambda²+36gamma²+64gamma²/K³]/pi².
The bound is deliberately conservative. Since the original zeroth tree wave is positive above threshold, the first elastic moment on any interval K>4 is strictly positive.

At K=10^198, use the actual rational couplings and pi>3. The resulting rational upper bound is strictly below gamma*10^-200. The conditional FULL lower bound in notes/split.md is greater than gamma/5. Therefore

J4_full/J4_first>2*10^199.

The conclusion is a required enhancement, not an assumption that such enhancement is impossible. A full absorptive calculation or matching model must determine its source, or one of the physical coefficient/analytic premises fails. The positive-pole control shows why it cannot be promoted to a whole-parent no-go.

For perspective, at COM energy at most10^99 the full original angular tree is bounded by4lambda K²+(3gamma/4)K³+8gamma<10^-202. Smallness of this tree value alone does not control the full moment. Nor does this first-cut estimate bound real quantum coefficients such as the mandatory b40.
