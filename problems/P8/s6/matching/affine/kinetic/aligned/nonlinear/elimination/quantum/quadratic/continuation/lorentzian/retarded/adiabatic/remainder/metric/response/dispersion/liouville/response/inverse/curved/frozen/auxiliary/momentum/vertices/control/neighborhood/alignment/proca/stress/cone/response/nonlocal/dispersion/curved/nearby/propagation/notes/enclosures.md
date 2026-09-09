# Exact whole-box cone bounds and the actual solution

The proof box is

    |u|<=10^-7,
    |N-(1+10^-6)|<=4*10^-7,
    |z|<=3*10^-5.

For each numerator and denominator polynomial, substitute
N=1+10^-6+n exactly, retain the exact constant term, and bound
every other monomial by its absolute rational coefficient times
the corresponding radius powers. Thus the polynomial lies in
center plus or minus the complete monomial sum. No Taylor
remainder or sampling assumption occurs: these are finite full
polynomials. If the resulting denominator interval avoids zero,
the minimum and maximum of the four endpoint quotients enclose
the rational function. Negative denominator intervals are allowed;
intervals meeting zero are rejected.

The exact rational endpoints are saved. Some weaker convenient
inequalities proved by those endpoints are

    5*10^-6 < c_clock^2-1 < 16*10^-6,
    -3 < PIV < -2,
    9/1000 < w^2 < 11/1000,
    -26/100 < M < -24/100.

The sharper saved speed endpoints also imply

    (1+2.5*10^-6)^2 < c_clock^2 < (1+8*10^-6)^2.

This last comparison is made with exact squares before taking
the positive characteristic speed; a speed-squared margin is not
mistaken for the same speed margin.

## Why the actual nearby bounce stays in this box

S6.88 proves, for this central lapse and on the entire interval,
|N-N0|<4*10^-7, |p|<=10^-6 and |log R|<10^-9.
The same actual coefficient b has a complex majorant below one
on |u|<=1/100, uniformly on the relevant lapse neighborhood.
It is odd in u, including its fixed-basepoint primitive. Cauchy's
coefficient estimate therefore gives

    |b| < rho/(1-rho^2), rho=10^-5.

This uses the complete odd analytic series, not only its first
time derivative. Since e<2,

    |z|=|e(p-b)| <2[10^-6+rho/(1-rho^2)]<3*10^-5.

For real u and N>1, hbg>=1 and
D-hbg=(N^2-1)(hbg-1)>=0, while N^2 hbg-D=N^2-1>0.
Hence 1<e<=sqrt(N)<1+10^-6 on the box. Also

    1-2*10^-9 <R<1+2*10^-9.

These scale inequalities follow, for example, from the exponential
series bound exp(x)<=1/(1-x) for 0<=x<1 and
exp(-x)>=1-x. They imply

    99/100 < N/(eR) <101/100.

All input radii come from the actual S6.88 solution. Eliminating
w^2 describes that solution because it satisfies its fixed-phase
constraint and M is nowhere zero. The positive square is also
explicitly enclosed, so the elimination is not continued onto an
unphysical negative-density branch.

## Positive principal forms and regular denominators

On the solution h=PIV/e<0 and m_N=eM!=0. In the triangular
characteristic basis of reduction.md, the exact bounds imply

    K_clock=-PIV/(4e^3 M^2)>7,
    K_matter=e^3/N>99/100.

Both physical characteristic values are positive, so both scalar
gradient eigenvalues in this basis are positive as well. The
unchanged tensor coefficients are F_T=G_T=2B4=e^-4>99/100.
The zero-field ordinary-Proca sector retains its unchanged
positive classical three-mode physical-metric principal system
from S6.81; no new quantum response is inferred.

These are positive principal forms on an actual short nearby
bounce, not infrared stability, a global nearby solution, nonlinear
stability or a uniform all-frequency canonical chart.
