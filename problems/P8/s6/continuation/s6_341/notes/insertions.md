# Literal TT line kernels and weighted bounds

For N scalar propagators with line shifts Qj and simplex weights xj,
barQ=sum xj Qj. Split line a into Qa and Qa+k, with weights xa-y and y.
The simplex Jacobian is1,0<=y<=xa; writing y=xa*gamma supplies the
essential factor xa*d_gamma. The exact second Symanzik polynomial changes by

 Delta_y-Delta=2y(barQ-Qa).k,

because k^2=0. This equals replacing Qa by Qa+gamma*k in the original
N-line denominator. The domain argument in parameters.md applies to every
gamma, including the endpoints and coincident-invariant limits.

The literal stress tensor from S295 is
p r^T+r p^T-eta(p.r-m^2), with r=p+k. Shift the loop variable by
barQ+y*k. Odd loop moments vanish; isotropic second moments contract
with trace(epsilon), which is zero in exact D. Transversality kills all
additional k terms. The numerator is therefore exactly

 2epsilon(Qa-barQ,Qa-barQ).

Before specializing D=4, the Feynman/radial coefficient is
Gamma(N+1-D/2)=Gamma(N-1+epsilon) and the denominator power is
N-1+epsilon. At D4 the triangle insertion has factor1 and power2;
the box insertion has factor2 and power3. The gamma integral is an
ordinary divided difference with its derivative limit at equality.
No longitudinal Ward reconstruction supplies this transverse numerator.

## Absolute bounds, with numerator-vector norm squared factored out

Use|xi' z'|<4, triangle measure<=2, box measure<=4t,
each light weight<=4 and each heavy weight<=2t. The full stress factor
is2; the box radial Gamma factor is2. With the moments of moments.md:

 triangle light: 8*4*2*64^2*32/n =8388608/n;
 triangle heavy: 8*2*2*64^2*500/n^2 =65536000/n^2;
 box light: 16*4*2*2*64^3*16/n^2 =1073741824/n^2;
 box heavy: 16*2*2*2*64^3*500/n^3 =16777216000/n^3.

These are bounds on weighted integrands after contour deformation; they
are not inferred by taking the modulus of an unweighted scalar result.

Since|xi|<=2x,|1-xi|<=2(1-x), the sum of absolute complex simplex
weights is<=4(1-t)+2t<=4. Original cyclic|Qj|<=12 gives
|Qa-barQ|<=60 and squared norm<10000. The canonical physical TT
projection is contractive in Euclidean Frobenius norm, so the same
bound controls both physical polarizations together.
Summing two light and one heavy line gives<10^12/n for a triangle;
two light and two heavy lines give<10^14/n^2 for a box.

These exclude the original couplings, common1/(16pi^2) loop factor and
1/sqrt(kappa). Overall Wick/vertex signs do not affect the absolute
bounds and are not being used to assemble a full radiative amplitude.
External emissions, outer heavy branches, labelled coefficients and
finite counterterms remain separate completion work.
