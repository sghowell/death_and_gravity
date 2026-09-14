# Full finite-momentum reference energy, with numerical constants

For E_y=(ydot^T K ydot+q y^T G y)/2 and the entire original
equation K yddot+(Kdot+3HK+gyro)ydot+(qG+lower)y=0,
antisymmetry removes only ydot^T gyro ydot. It does not remove
gyro from the canonical momentum or chart maps. The two diagonal
energy numerators are-(Kdot+6HK) and Gdot-2HG; the remaining work
is-ydot^T lower y. Every lower term is retained.

On each of the three source-pinned chart intervals subdivide into64
closed rational cells. For the full profile-dependent principal K0,G,
evaluate the first principal minor and determinant of
20Q+V,20Q-V,Q-I/2 and32I-Q for Q=K0,G and their respective full
energy numerators V. These are3072 outward rational conditions,
including both nonzero original profiles and first jets.
Every lower endpoint is greater than1/100. Sylvester's criterion
gives K0,G in[1/2,32] and absolute principal root-energy rate10.
Sampled eigenvalues are separate diagnostics, not this proof.

Only the central kinetic matrix requires a finite-q correction.
Put z=1/q. Each entry of K-K0 and its ENTIRE time derivative is
normalized by E^4 Delta^2/z and reconstructed exactly as a rational
polynomial. The old valid coefficient/jet box, including qdot=-2Hq,
bounds every normalized entry below10^40. E is bounded away from
zero and Delta>=1/8. Hence each matrix difference has operator
norm at most eta=2*10^40/q_min, q_min=10^128/4.
All eight polynomials are checked; this is a bound at each finite
P>=10^64, not an asymptotic inference.

The source-pinned original lower matrix entry bound is10^18,
including the full profiles, lower potentials and coefficient jets.
Its operator norm is at most2*10^18. With the safe kinetic/gradient
floor1/4 its work contribution is bounded by
2*10^18/[(1/4)sqrt(q_min)], less than1/100.
For the diagonal numerator, the LMI perturbation budget
(20+1+12)eta includes20 times the K change, the Kdot change
and6|H| times the K change, with |H|<=2.
It too is less than1/100. The positive principal margin survives:
the full finite-q K,G lie in[1/4,64].
Converting the small work/numerator errors through this floor gives
the deliberately loose absolute full root-energy rate12.

The clean canonical momentum is Pi=a^3(K ydot+A_skew y);
both symmetric boundary terms are kept in the separate chart maps.
For x=(P y,Pi), reconstruct the full symmetric energy metric Q from
ydot=K^-1(x_p/a^3-A_skew x_q/P) and
2E_y=ydot^T K ydot+x_q^T G x_q/a^2=x^T Q x.
Both full chart coefficient bindings and all mixed blocks are stored.

Use1<=a<=2, K,G in[1/4,64] and ||A_skew||<=2*10^18.
Young's inequality gives a lower eigenvalue at least
min(1/(128*64),1/16-(2*10^18/P)^2/64)>10^-4.
The upper bound is at most8+8(2*10^18/P)^2+64<100.
No commutativity or simultaneous diagonalization of K and G is assumed.

Finally r=(sqrt(P)y,Pi/sqrt(P)) is canonical, whereas x=sqrt(P)r.
Thus E_y=P r^T Q r/2 and Q in[10^-4,100].
Using x as if its CCR were Omega would lose a required momentum
factor in every subsequent covariance and physical point bound.
