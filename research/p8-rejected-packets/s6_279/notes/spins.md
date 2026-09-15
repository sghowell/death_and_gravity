# All-even-spin proof, beyond the finite exact examples

For0<beta<1 put w=1/beta and C=s(1-beta^2)^2/(4kappa).
The normalized scalar and helicity-difference4 basis functions are
P_l(x) and sqrt((l-4)!/(l+4)!)*(1-x^2)^2 P_l''''(x).
Each has squared integral2/(2l+1). The latter exists only for l>=4.

The real-sheet Legendre function Q_l(w) is
one half of integral_-1^1 P_l(x)/(w-x) dx, w>1.
The same-helicity coefficient follows by partial fractions:

    f_l0=C(2l+1)Q_l(w)/beta.

For the opposite-helicity coefficient, four integrations by parts move
the derivatives from P_l onto (1-x^2)^4/(w^2-x^2).
That numerator and its first three derivatives vanish at both endpoints.
Polynomial division leaves a degree6 quotient and remainder(1-w^2)^4.
The fourth derivative of the quotient has degree2, hence is orthogonal
to EVERY l>=4, not merely the checked examples.
Differentiating the Q integral then gives

    f_l4=C(2l+1)sqrt((l-4)!/(l+4)!)(w^2-1)^2 Q_l''''(w)/beta.

For l<4 set f_l4=0. Odd l vanish by parity. The code verifies the Q
integral, both norms, remainder and full integrated coefficient through
l12, but the preceding argument establishes the formula for all spins.

Sewing one fixed physical helicity uses the spin-weighted harmonic
addition theorem. Orthogonality over the cut direction gives
P_l(z)/(2l+1). There are two same and two opposite helicities; combining
this with the angular32pi normalization yields

    Im A(s,z)=sum_even_l
      [f_l0^2+f_l4^2] P_l(z)/(16pi(2l+1)).

All weights are nonnegative because Q_l and its derivatives are real on
w>1. The explicit amplitudes are analytic smooth spin sections on the
compact sphere for beta<1: their denominator is strictly positive, and
the opposite-helicity numerator has the required fourth-order polar
zero. Harmonic coefficients therefore decay faster than any fixed
power (in fact exponentially at fixed beta). The sewn series converges
absolutely. Parseval at z1 agrees with the independent rational forward
integral, and |P_l(z)|<=1 controls every physical angle.

No finite list of verified polynomials is presented as an infinite-spin
computation. This argument proves positivity of this perturbative
physical cut, not UV admissibility of an arbitrary smearing functional,
the exact finite-detector amplitude or a Regge bound.
