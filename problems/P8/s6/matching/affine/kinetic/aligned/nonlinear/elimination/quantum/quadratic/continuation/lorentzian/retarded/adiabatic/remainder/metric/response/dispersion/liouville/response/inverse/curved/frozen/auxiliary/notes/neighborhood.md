# A convergent auxiliary branch with explicit constants

All derivatives in this note act on the actual Hcal of
hamiltonian.md, not a finite lapse jet substituted for it.
The classical tree plus margin and the exact spatial
canonical-jet invariants are fixed.

## Uniform Cauchy estimates and a closed-disc map

Fix real u in I. The outer holomorphic bound is
|Hcal|<M=10000 for |N-1|<=1/100 and max|Y_i|<=1.
On the inner set |N-1|<=1/200, max|Y_i|<=1/2,
Cauchy's formula with remaining radii 1/200 and 1/2 gives

|Hcal_NYi| <= 4000000,
|Hcal_NNN| <= 480000000000,
|Hcal_NNYi| <= 1600000000.

The factorial 2 is retained in the last bound, and the
factorial 6 in the middle bound. Write F=Hcal_N and
d0=Hcal_NN(1,0)=-2J_e, with |d0^-1|<20.

Take RN=10^-14 and RY=10^-24. For |n|<=RN and
max|Y_i|<=RY, the straight line from (1,0) to (1+n,Y)
remains in the inner polydisc. Therefore

|F_N(1+n,Y)-d0|
 <= 480000000000*RN+9*1600000000*RY
 = 3000000000009/625000000000000,

|F(1,Y)| <= 9*4000000*RY.

The second line uses the exact background identity
F(1,0)=0. The bounds are for complex inputs too.

For each Y, use the Newton map with fixed background
derivative,

Phi_Y(n)=n-d0^-1*F(1+n,Y).

Its derivative has modulus at most
3000000000009/31250000000000 < 1/10.
Its value at n=0 has modulus at most
20*9*4000000*RY.
The sum of this center step and the contraction bound
times RN is strictly smaller than RN, as checked exactly.
Thus Phi maps the closed complex n disc into itself and
is a strict contraction there.

Banach's theorem gives one fixed point, and every zero
of F in this disc is that fixed point. Iteration from zero
is holomorphic in Y and converges uniformly on the closed
polydisc; the root is holomorphic on its interior.
The strict margins and holomorphic original coefficients
also extend the construction to an open neighborhood of
the closed polydisc. This justifies Cauchy estimates on
its boundary without an endpoint assumption.

For an actual input radius r=max|Y_i|<=RY, the same
bounds yield

|n| <= (10/9)*20*9*4000000*r = 8*10^8*r.

At r=RY this is 8*10^-16, below RN/10 and below 10^-15.
The secondary derivative remains bounded away from zero:

|Hcal_NN| >= 1/20
 -3000000000009/625000000000000 > 1/25.

For real u,Y, all coefficient branches respect complex
conjugation. Uniqueness then makes the root real. Its
N is positive and X=N^-2 remains strictly inside
[9/10,11/10]. The implicit branch varies smoothly in u
by the original smooth coefficients and nonzero pivot.

## Temporal reconstruction and the full auxiliary block

Let G=K-3H/N. From hamiltonian.md, |K|<27 and
|3H/N|<18 on the outer N/Y polydisc, so |G|<45.
The exact background reconstruction gives G(1,0)=0.
The same inner Cauchy radii give

|G| <= 9000*|n|+9*90*r.

For |n|<=RN, the exact N-to-s map gives
|s^2-1|<=3|n|, hence |delta|<=3|n| for real u.
The original Q estimate and |s|<2 imply

T=delta*G+(3/2)*s*Q-gamma_t*j/U,

|T| <= 3|n|*(9000|n|+810r)+108|n|^2+4r.

Substitute |n|<=8*10^8*r and r<=10^-24. The exact
rational coefficient is at most
502168640000243/125000000000000 < 5.
Thus |T|<5r and in particular |T|<5*10^-24.
The source's zero first variation is used here; it has
not been discarded nonlinearly.

S6.45 proves gamma_eff>9/10 throughout the original
real tube, so the temporal auxiliary Hessian before its
elimination is nonzero. The derivative Hcal_NN just
bounded is the Schur complement after temporal
elimination. Their product gives the determinant of the
full joint auxiliary Hessian, which is therefore nonzero
on this explicit real canonical-jet domain.

The remaining ten velocity pivots and the spatial
first-class identities are the full S6.45 ones, unchanged
by the scalar margin. Consequently its local
seven-degree-of-freedom count persists here. This is
local constraint solvability/count, not hyperbolic
evolution or nonlinear energy stability.

## Invariant Taylor coefficients and a true remainder

Substitute the holomorphic root into the full Hcal.
The result Hred(Y) stays below 10000 on the closed
nine-invariant polydisc of radius RY. Multivariate Cauchy
estimates give, for ordinary monomial coefficients,

abs([Y^alpha] Hred) <= 10000*RY^(-|alpha|).

There are binomial(k+8,8) monomials of total degree k.
Therefore the total coefficient bounds at degrees three
and four are 165*10000*RY^-3 and 495*10000*RY^-4.
This is an invariant Hamiltonian bound, not a
momentum-reduced scalar/tensor/vector vertex calculation.

For max|Y_i|<=theta*RY, theta=1/1000, the tail after
degree four is bounded by

10000*sum_(k>=5) binomial(k+8,8)*theta^k.

The ratio of consecutive summands is
theta*(k+9)/(k+1), at most 7/3000 for k>=5.
Its first summand is 10000*1287*10^-15, giving

abs(remainder) <= 3861/299300000000 < 10^-7

on max|Y_i|<=10^-27. This is a continuous-domain
remainder of the entire actual implicit Hamiltonian,
not comparison at selected sample points. Its physical
normal-density scale is restored by M^2/tau^2 in the
same fixed hat-volume convention.
