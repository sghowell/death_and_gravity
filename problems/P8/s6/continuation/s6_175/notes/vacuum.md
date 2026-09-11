# Complete analytic source: a full-function small-field bound

Use the original real-spacetime test class, not a Wick-rotated claim:
Phi is real Schwartz on R^4, its Fourier support lies in the Euclidean
four-momentum unit ball, and its Fourier L1 norm is at most one in the
fixed inverse-transform convention. Put U=||Phi||2 and u=Phi/sqrt(kappa).
Every component of every derivative jet of Phi has supremum at most one
and L2 norm at most U, by Fourier inversion and Plancherel.
The word Euclidean in the support ball specifies a momentum norm only.

Let even n>=4, kappa>=64n and rho=1/(16n). In the actual target n=1024,
kappa=10^800. Then |X|<=4/kappa<=rho. For the COMPLETE functions

    T=X^n/[X^n+(1-X)^n],  w=n X^2 exp(-n X^2),
    B=T+(1-T)w,  h=(1+u^2)^3,  R=1+B(X-1)/h,

one has 0<=T,w<=1. Since |X|<=rho<1/2,

    T <= (2|X|)^n,
    T/X^2 <= 4(2rho)^(n-2) <= 1,
    |T'| <= 8n|X|(2rho)^(n-2) <= n|X|,
    |w'| <= 2n|X|.

The w' inequality uses 0<=nX^2<1. The displayed limits at X=0 are
removable. Therefore B<=(n+1)X^2 and |B'|<=3n|X|.
Keeping the factor |1-X| is important: it is bounded by 1+rho, not one.
The elementary positive margins checked in the certificate imply

    |R-1| <= 2n X^2/h,       |r| <= 2n/h,
    |R_X| <= 4n|X|/h,       R >= 1-2n rho^2 > 1/2,
    r=(R-1)/X^2,            R_u=-6u(R-1)/(1+u^2).

These estimates hold for positive, negative and zero X. At n=1024
the interval [-rho,rho] lies strictly inside the parent analytic domain.
With H=4u/(1+u^2), the source is evaluated in its regular form

    S_mu=u_mu [(R-1)(-3H+3R_u/(4R))
               +X r Box u +(-r+3X r R_X/(2R)) Z],
    Z=u^mu u^nu u_;mu nu.

There is no division by X in this form. The jet bounds give

    |Box u|<=4/sqrt(kappa), |Z|<=16/kappa^(3/2),
    |R-1|<=32n/kappa^2, |R_X|<=16n/kappa,
    |R_u|<=192n/kappa^(5/2), |H|<=4/sqrt(kappa).

Place the outer u_mu in L2 and the remaining factors in L-infinity.
For EACH component the resulting coefficient multiplying U is

    64n/kappa^2 +384n/kappa^3
      +6144n^2/kappa^4 +9216n^2/kappa^5
      <=128n/kappa^2.

To check the last inequality, divide by n/kappa^2, replace kappa by
64n in the positive inverse powers, then n by4. The result is
269833/4096<128. Taking the ordinary positive four-component L2 norm
gives ||S||2<=256n kappa^-2 U. No infinite spacetime volume is used,
and no Taylor truncation of R, its derivatives or the source is made.

For an INDEPENDENT real Euclidean functional, define the continued
scalar contractions X_L=-X_E, Box_L=-Box_E, Z_L=Z_E and its real
Euclidean one-form source. The same absolute estimates apply.
O_E=I+zeta(p_E^2 I-p p^T) is positive and

    0 <= I-O_E^-1 <= I.

At p=0 this follows by continuity, without a singular longitudinal
projector. The full source Schur functional consequently satisfies

    0 <= Delta S_E <= kappa ||S_E||2^2/2
      <=32768 n^2 kappa^-3 U^2 <10^-2389 U^2.

This is a full-source Euclidean observable. It is NOT an analytic
continuation theorem for the original Lorentzian real test class,
a real-time matching estimate, or an interacting parent quantum bound.
