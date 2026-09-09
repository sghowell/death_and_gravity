# Actual complex-time bounds and uniform normal-form constants

Write T=10^-7. S6.88 gives the actual analytic solution on the
complex disc |u|<=T, not only the real interval. Its phase image,
lapse-root displacement and scale reconstruction fit S6.89's
complex coefficient box. The odd analytic bound for b and
|e|<2 give |z|<3*10^-5 on this same solution.

For a full polynomial with exact center a and complete monomial
deviation bound d, its complex modulus lies between
max(0,abs(a)-d) and abs(a)+d. Applying this to numerator and
denominator separately gives genuine complex annuli; a real
interval inequality alone is not treated as complex ordering.

The exact saved annuli imply

    5e-6<abs(c_clock^2-1)<16e-6, 2<abs(PIV)<3,
    .009<abs(w^2)<.011, .24<abs(M)<.26, abs(alpha)<1.

With D=1+N^2[(1+u^2)^3-1],

    e^4-1=(N^2-1)/D.

The complete binomial-series majorant gives |e-1|<4e-6.
Likewise the exponential-series majorant for the actual log R
gives |R-1|<2e-9. Both e and R therefore lie in the .99-to-1.01
modulus annulus, and |N-1|<1.4e-6.

The clock root c=sqrt(1+excess) is the branch near1, with
|c-1|<1/2. Thus 1/2<|omega_m|<2 and
1/4<|omega_j|<4. The clock/matter same-sign gap has modulus

    |omega_m| |excess|/|c+1| >(.5)(5e-6)/3 >10^-7.

Opposite-sign gaps use Re(c)>0 and |c+1|>1.
The gaps 2omega_j are also separated. All four signed frequency
gaps are therefore uniformly greater than 10^-7 in modulus.
The analytic nonzero kinetic factors have square roots on the
simply connected actual time disc, fixed by their real positive
central values; no root branch is changed along the solution.

## Explicit coefficient and basis majorants

The literal identities and annuli give the weak bounds

    1<|h|<6, .1<|r_N|<1, |ell|<1/8,
    |a|<3, |m|=|r|<8, |g|<2,
    |alpha|,|beta|,|Hhat|<1,
    .5<|kappa_c|<200, .01<|kappa_m|<100.

The program expands each finite Laurent coefficient and sums its
absolute monomial bound, using lower moduli for negative powers.
Every entry of A0,A1,A2,E0,E1,C0 is below100.
Including |R|<2 and the volume terms gives ||Lj||_infinity<10^4.

The explicit mode matrices are Laurent polynomials in the positive
kinetic/frequency square roots, analytically continued on the disc.
Their complete entry bounds give

    ||S||_infinity <=85sqrt(2)/2 <1000,
    ||S^-1||_infinity <=89sqrt(2)/4 <1000.

For complex time these are modulus majorants of the same analytic
formulas, not complex-conjugate differentiations.

## Time Cauchy bounds on an inner real interval

Work on I=[-T/2,T/2]. Every point has complex-time Cauchy radius
T/2 within the actual solution disc. Therefore

    ||S'||,||(S^-1)'|| <=2e10,
    ||S''|| <=8e17, ||L0'||<=2e11.

The second derivative retains its factorial. These are derivatives
along the actual solution, not loose bounds formed by multiplying
independent phase derivatives. Product rules now give

    ||B0||<3e13, ||B0'||<2e21, ||Lambda'||<=8e7.

For k>=1,

    W=S^-1(L1+L2/k+L3/k^2)S, ||W||<3e10.

Set Z_ij=i B0_ij/(lambda_i-lambda_j), Z_ii=0.
The exact gap equation and its derivative, including derivatives
of the gaps, give ||Z||<1e21 and ||Z'||<1e36.
All displayed rounded constants are weaker than saved exact
rational triangle bounds.
