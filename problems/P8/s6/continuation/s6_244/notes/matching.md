# Complete general-dimensional scalar matching

## Universal one-oscillator variation

For arbitrary positive K and omega, let lambda=omega'/omega, s=(lambda-K'/K)/2 and t=s'-lambda s. The universal S193 identity is used ONLY in its single-oscillator specialization:

L0=-omega/2, L2=s²/(4omega), L4=(t²+s⁴)/(16omega³).

Its complete compact Euler variations equal the physical diagonal covariance-current coefficients through adiabatic order4. Here this is independently checked again using raw scalar omega0..4 and eta0..4, eta=log K: the literal Euler derivatives of each action are compared to QQ and PP from the full scalar Riccati recursion. All six independent variations vanish exactly. Neither a vector multiplicity nor a vector finite counterterm is inherited.

## Arbitrary homogeneous spatial metric

Work in a Fermi orthonormal frame. Let Kgeom be the symmetric extrinsic curvature, A=D_t Kgeom, theta=tr Kgeom, X=tr A. No commutation of Kgeom and A is assumed. Let u be the physical unit momentum direction, z=q²/(n+q²), and k_u=u^t Kgeom u. Then

lambda=-z k_u,
(k_u)'=u^t A u-2u^t Kgeom² u+2k_u²,
s=(theta-z k_u)/2,
t=[X+z u^t(-A+2Kgeom²+theta Kgeom)u-3z² k_u²]/2.

The complete physical spatial curvature is zero, while the spacetime invariants are

R=2X+theta²+k2,
Ric²=(X+k2)²+A2+2theta AK+theta²k2,
Riem²=4A2+8AK2+2k4+2k2²,

where kj=tr Kgeom^j, A2=tr A², AK=tr(A Kgeom), AK2=tr(A Kgeom²). Literal Christoffel/Riemann contractions independently test these formulas at dimensions2,3,4,5 with two genuinely noncommuting Kgeom,A fixtures per dimension.

For general spatial dimension d, the full angular moment of N quadratic forms has Wick denominator d(d+2)...(d+2N-2). Its radial z^N moment relative to the base omega^-2alpha integral contributes (d/2)_N/(alpha)_N. The product has denominator (2alpha)(2alpha+2)..., independent of d, without setting the metric dimension to3. The complete trace invariants themselves retain their dimension dependence.

At adiabatic order2, alpha=1/2. The full average of (theta-z k_u)² is2(k2-theta²)/3. Together with the weighted boundary R-2v^-1(v theta)'=k2-theta² this yields the whole Einstein density.

At order4, alpha=3/2. Expanding EVERY term of t²+s⁴ and applying the full moments gives

5theta⁴/756+theta²X/21+16theta²k2/315-theta AK/105
-34theta k3/945+X²/10+2Xk2/21+A2/30
-2AK2/105+7k2²/180+k4/45.

Its difference from P=R²/36+(Riem²-Ric²)/90 is exactly

-4/189*(theta⁴+3theta²X)
+2/315*(theta²k2+Xk2+2theta AK)
-34/945*(theta k3+3AK2).

These are v^-1 derivatives of v theta³, v theta k2, v k3. They are full weighted boundary terms at GENERAL dimension and vanish under compact metric variations before taking the dimensional limit. They are not discarded only after d=3.

Cartesian polynomial angular integration plus the exact radial ratio independently tests the whole fourth-order polynomial in all eight noncommuting fixtures. An arbitrary-d isotropic reduction also agrees exactly.

## Entire MSbar normalization and finite action

Set d=3-2epsilon, ell=log n and F(epsilon)=exp[(gamma_E-ell)epsilon]Gamma(1+epsilon). Factoring the COMMON full dimensionally continued volume and invariants, and normalizing the action by64pi², the complete radial coefficients multiplied by epsilon are

2F/[(epsilon-1)(epsilon-2)] for n²,
F/[3(epsilon-1)] for nR,
F for P.

Their pole values are(1,-1/3,1); their finite values are(3/2-ell,(ell-1)/3,-ell). Thus the full pole action is n²-nR/3+P, with P=2a2 for the minimal scalar.

Crucially, the counteraction retains the SAME dimensionally continued v,R,Ric²,Riem² before variation and the limit. If a full invariant is H(d)=H0-2epsilon H_d+..., the finite part of

[(coefficient(epsilon)-pole)/epsilon] H(d)

is precisely finite_coefficient*H0: the potentially nonzero evanescent derivative cancels between the full action and full counteraction. The code differentiates the complete expression, not isolated physical components. Weighted total derivatives from the preceding section are already zero for compact variations at general d.

The resulting complete physical finite density is EXACTLY the existing S240 prescription:

[(3/2-ell)n²+(ell-1)nR/3-ell*(R²/36+(Riem²-Ric²)/90)]/(64pi²).

No renormalization parameter is refitted. The formula holds for the entire homogeneous metric family and hence all compact first/second variations considered here. Its local current is added ONCE to the convergent actual-state-minus-adiabatic integral.

S243 concerns a spatial transfer DIFFERENCE, which is zero at this anchor. Its finite difference is not an extra homogeneous term. Likewise S241's finite heat Hessian is the variation of this same action, not another action to add.
