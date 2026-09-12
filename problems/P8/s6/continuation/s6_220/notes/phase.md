# Regular phase comparison and its limits

Add the full fixed retuning quadratic to the literal action. The normalized momenta are
\[
 p_v=-6\dot v+6\Theta n-3\ell\sigma+2b,\qquad
 p_\sigma=\dot\sigma+wn.
\]
Perform the Legendre transform before eliminating the two auxiliary variables. Its b derivative is(p_v-2b)/3, so b=pv/2. After that substitution its n derivative is-2Jc n+Lc. Thus n=Lc/(2Jc), with the Lc and Hamiltonian in FORMULATION.md. All identities are exact at Theta=0.

The normalization by a^3 makes the equations
\[
 Z'=J_{\rm can}\operatorname{Hess}{\cal H}_c\,Z
       -3H\operatorname{diag}(0,0,1,1)Z.
\]
The weighted momentum term cannot be omitted. Independent direct linear saddle solves verify both momenta and constraints at positive, negative and exactly zero Theta, with very small and large q. They do not use the displayed lapse or shift formulas to solve the saddle.

This coefficient-sector Hamiltonian is not automatically coercive. At the bounce, choose q=1, v=ps=0, sigma=1,pv=20. The fixed retuning terms and lapse square vanish for this direction, and the Hamiltonian is-203/400. This is a control against a positive-energy overclaim. It is not a conclusion about the full time-dependent physical scalar spectrum; time-dependent canonical changes may alter an instantaneous Hamiltonian.

## Explicit coefficient and compact-momentum bounds

The clock formula gives |Theta|<=|H|<=2. Also ell<=1/10, |w|<=1/10, |E|<=1. With |A|<epsilon<1/100, |Tc|<(5/2)epsilon<1/100 and Jc>1/100, the coefficient row of Lc obeys
\[
 \|L_c\|_{\rm row}\le2q+273/100<3(1+q),\quad
 1/(2J_c)<50.
\]
The base symmetric Hamiltonian Hessian has row bounds
\[
 2q+9/100,\quad q+13/200,\quad1/20,\quad1,
\]
each below3(1+q). Its operator norm is bounded by its maximum absolute row sum. The rank-one lapse Hessian is l l^T/(2Jc), bounded by450(1+q)^2. Since q<=|P|^2 and the weighted damping norm is<=6,
\[
 \|K(t,P)\|_2<1000(1+|P|^2)^2.
\]
The exact margin before replacing q by |P|^2 is541+1097q+550q^2, strictly positive.

At each0<|P|<=Lambda the smooth matrix IVP has its unique fundamental solution. Volterra iteration, bounded by the exponential series, yields
\[
 \|U(t,s;P)\|_2\le
 \exp[1000(1+\Lambda^2)^2|t-s|].
\]
For zero comparison initial data, the retarded variation is integral_s^t U(t,r)f(r)dr. The unit interval and Schur/Young bound give L2-time norm no larger than the same slab exponential. This extends as a bounded Fourier multiplier on the compact-momentum R3 subspace. The value at the origin is immaterial there; Lambda=0 gives the trivial R3 subspace, not a new homogeneous constraint result.

This is a classical coefficient-sector propagator. It neither resets the quantum state nor solves the full Gaussian lapse constraint. Its exponential grows as Lambda^4; it is not a uniform continuum estimate or a choice of physical EFT cutoff.

## Actual coefficient-sector metric map and derivative loss

With the actual smooth fixed QG1 coefficients, write the linear map
\[
 (n,\zeta,b)=L_0(t)Z+|P|^2L_2(t)Z,
 \quad n=L_c/(2J_c),\quad\zeta=v+\delta n,\quad b=p_v/2.
\]
The implementation constructs both matrices from this exact map, including a^-2 in L2. The preceding bounds give
\[
 |n|\le200\lambda|Z|,\quad|\zeta|\le101\lambda|Z|,
 \quad|b|\le\tfrac12\lambda|Z|,
\]
so V03[sD]<=400 V05[ZD].

For source jets, define the explicit finite number
\[
 C_{13}=\sum_{r=0}^{13}\sum_{j=0}^r{r\choose j}
 \sup_{t\in I}\left[
  \sum_{\alpha,\beta}|(L_0^{(j)})_{\alpha\beta}|
 +\sum_{\alpha,\beta}|(L_2^{(j)})_{\alpha\beta}|\right].
\]
Smoothness of the fixed reference coefficients and Jc>1/100 make this finite. Each entrywise matrix norm bounds its operator norm. Apply the product rule to each time jet and sum the resulting finite block-operator norms. This yields U138[sG]<=C13 U13,10[ZG]. It does not infer thirteen small stress-jet bounds from the established five-jet budget.

Combining with the infrared scalar metric form gives
\[
 |R_{\rm phase}(D,G)|<1e120 C_{13}V_{05}[D]U_{13,10}[G].
\]
The source has thirteen time and ten spatial derivatives; the detector five spatial derivatives. Both chart contacts are retained. This is the response along a fixed linear comparison map, not the Hessian of a fully quantum-eliminated nonlinear reduced action.

The full Gaussian lapse equation has nonlocal response. A genuinely coupled inverse/Schur argument is still required. No canonical1/kappa factor is transferred to unspecified reduced norms or source variance, and these derivative-losing bounds do not provide a direct same-space Neumann contraction.
