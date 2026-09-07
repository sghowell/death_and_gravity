# A quantitative full physical loading Gramian

All vector and matrix norms here are Euclidean unless an infinity norm is
explicitly named. The source Hilbert norm is L2 in u, not in the rescaled
coordinate. This proof uses the pinned S6.23 and S6.21 operator, source,
wave and parameter-transfer theorems; it changes no frozen ancestor.

## 1. Actual source and a regular scaled physical system

Set ell=1/100, a=-2ell, b=-ell and x=(u-a)/ell. Only on this punctured
interval may we set delta=0 in the coefficients. This is not a new
nonsingular action at the bounce. Write d=1+u^2 and K in [1,4]. The
literal unit-polarization equations are

    (Kg g_u)' + Gg K g + U(g-f) = a_scale^3 sigma/4,
    (Kf f_u)' + Gf K f + U(f-g) = 0,

where a_scale=d^2, Kg=d^6/8, Kf=1/(2d^6),
Gg=d^2/8, Gf=1/(2d^2), and

    U = 4(1-u^2)/[d^8(d^4-1)] > 0.

Thus the physical derivative-state source is exactly (0,2,0,0)sigma.
The actual canonical-to-physical map and its derivative reproduce all
16 generator entries and all four source entries in core.checks().
In particular an arbitrary relative source has not been substituted.

For X=(g,ell*g_u,f,ell*f_u), the x generator N is

    [       0,          1,          0,          0 ]
    [ -ell^2 K/d^4-Ug, -p,         Ug,          0 ]
    [       0,          0,          0,          1 ]
    [      Uf,          0, -ell^2 K d^4-Uf,     p ],

with Ug=ell^2 U/Kg, Uf=ell^2 U/Kf, p=12ell*u/d. Its input is
2ell^2 e2 sigma(a+ell*x). Every denominator is positive or bounded away
from zero on this closed interval. Factoring
d^4-1=u^2(4+6u^2+4u^4+u^6) also avoids interval cancellation.

Let G(0)=0 solve G_x=NG+GN^T+e2 e2^T. The physical controllability
Gramian for the L2(du) source is 4ell^3 G(1), since du=ell dx.
Using 4ell^4 here would give the wrong source cost.

## 2. Validated midpoint integration, not a sampled verdict

gramian.py integrates this exact rational-coefficient matrix ODE at K=5/2.
The primary replay uses 64 adjacent closed steps of length h=1/64,
Taylor order18, and 256-bit Arb ball arithmetic. Matrix products, series
reciprocals, roots in the final map and comparisons are outward-enclosed.
No floating-point trajectory, fitted error or accepted ODE-solver status
is used as certificate evidence.

Here is the complete remainder argument. At a step's left endpoint,
the current ball matrix encloses the exact G0. Evaluate the Taylor
arithmetic for N on the entire closed step, with a formal unit first
derivative. Its coefficient N_j encloses N^(j)(xi)/j! for every xi in
that step. This works because each denominator excludes zero there.
The norm cap

    lambda >= ||N||_infinity + ||N^T||_infinity

holds continuously on the whole step, and the replay checks h*lambda<1.
The integral equation then bounds the exact state throughout the step by

    ||G(xi)||_infinity <= (||G0||_infinity+h)/(1-h*lambda).

Every entry of a whole-step state is put inside the resulting centered
ball. Differentiating the ODE gives the exact normalized-derivative
recurrence

    G_(n+1) = [ sum_(j=0)^n (N_j G_(n-j)+G_(n-j) N_j^T)
                  + 1_(n=0) e2 e2^T ]/(n+1).

Applying this recurrence to the whole-step N jets and whole-step state
encloses G^(19)(xi)/19!. Entrywise Taylor's theorem bounds the remainder
by its absolute enclosure times h^19. The endpoint Taylor polynomial
uses the separate left-endpoint jets and the previously enclosed G0.
Adding the remainder balls closes the induction on all 64 steps.
This proof bounds the derivative of the actual solution, not merely a
truncated-series residual or a pointwise coefficient scan.

The primary replay records a step lambda below26 and a largest local
entry remainder below 1/10^27 (the checked general cap is 1/10^18).
Exact rational outward rounding of the final report matrices is applied
only after the calculation. The runtime restores Arb precision and
series-cap settings even if a validation check raises an exception.

For the final chart set r=-u, mu=sqrt(39)/2, and

    Z = (l, ell*l_u, Q/sqrt(r), (u Q_u-Q/2)/(mu sqrt(r))).

The endpoint map X->Z retains derivatives of every coefficient in
l=fs(w1*g+w2*f), Q=fr(f-g). Here
fs=sqrt(2d^12+8)/(2sqrt(2)d^3),
fr=sqrt(2)d^3/sqrt(2d^12+8), w1=2d^12/(2d^12+8), w2=8/(2d^12+8).
The minus-side orientation in the fourth row is essential.

After this exact map, interval LDL verifies all four pivots of

    W_Z(5/2) - (36/10^12) I

are strictly positive. The enclosed exact matrix is symmetric by its
ODE/Gramian definition; no symmetry of independent ball radii is needed.
Therefore its least source-map singular value is greater than 6/10^6.
The certificate records rational intervals for all four pivots and the
full matrices, not rounded eigenvalues.

## 3. A continuous momentum theorem

A midpoint Gramian alone is not a statement on [1,4]. To extend it use
the actual Z chart, not a naive norm bound on the stiff physical matrix.
In x time its exact generator is

    [ 0, 1, 0, 0 ]
    [ -ell^2 A, 0, -ell^2 sqrt(r)(C+dc/2), -ell^2 sqrt(r)dc mu ]
    [ 0, 0, 0, -ell mu/r ]
    [ ell sqrt(r) E/mu, -r^(3/2) fc/mu,
          ell mu/r+ell r b_rem/mu, 0 ],

where b_rem=B-10/r^2 is the punctured bounded remainder. The chart is
not defined by assuming that this remainder is jointly continuous at
the excluded bounce origin. A generic chain-rule identity verifies
every displayed entry and the input

    Bx = (0, ell^2 jL, 0, -ell sqrt(r) jH/mu).

The S6.21 bounds |A|,|C|<28, |dc|,|fc|<14,
|E|<=160r^2, |b_rem|<200 and |jL|,|jH|<1 hold here.
Separate the exact heavy skew rotation and the light entry N12=1.
The latter has symmetric-part norm1/2. Using
r<=1/50, sqrt(r)<1/7 and 3<mu<13/4, the entry-sum norm of all the
remaining entries is strictly below1/20. Consequently the logarithmic
energy exponent is bounded by kappa=11/20, and on any subinterval

    ||Phi_K(t,s)|| <= exp(kappa*(t-s)),       0<=s<=t<=1.

The elementary factorial bound
exp(q)<=1+q+q^2/[2(1-q/3)] for 0<=q<=2/3 proves exp(kappa)<2.

The exact momentum dependence is affine. Since
d^4<101/100 and sqrt(w1w2)<=1/2,
the absolute K derivatives of A and B are at most d^4<101/100.
The shared derivative of C and E is bounded by
(d^4-d^-4)/2 < d^4-1 < 1/100. The other listed entries and Bx
are K-independent. The sum of all K-derivative entry bounds is
strictly below1/5000. Also

    ||Bx||^2 < ell^4+ell^2*(1/50)/9 < (1/2000)^2.

Duhamel's formula, with the duration-dependent exponential bound, yields

    ||Phi_K(1,t)-Phi_(5/2)(1,t)||
        < 2 |K-5/2| (1-t)/5000.

The two propagator factors multiply to exp(kappa*(1-t))<2, not4;
this use of their durations is part of the proof. Integrating the
source kernel and converting its L2 norm from dx to du gives

    ||C_Z,K-C_Z,5/2||
      <= 2*(1/2000)*(1/5000)*(3/2)/sqrt(ell)
       = 3/10^6

uniformly for the whole closed momentum interval. Dropping the further
factor 1/sqrt(3) from the (1-t) integral only enlarges this bound.
The least singular value perturbation inequality now proves

    W_Z(K) >= (9/10^12) I                     for every 1<=K<=4.

All constant comparisons are exact rationals in core.py and are
separately reconstructed in independent.py.

## 4. Return to the actual target frame

Let C_K sigma denote the full four incoming coefficients in the
S6.21 Frobenius-normalized wave frame at u=-ell. Its physical definition is

    C_K sigma = Psi_-(ell)^-1 W_-(ell)^-1
                   (l,l_u,Q,Q_u)_sigma(-ell).

It is NOT raw Y and does not identify Y with the light inclusion Lv.
Writing D_ell=diag(1,ell,1,1), the above source map obeys

    C_Z,K = A_K C_K,       A_K=D_ell Psi_-(ell).

The pinned ||Psi_-(ell)||<2 implies ||A_K||<2. All maps are invertible,
so W_Z=A_K W_C A_K^T and

    W_C(K) >= [9/(4*10^12)] I.

The pinned source kernel bound ||C_K sigma||<=8||sigma||_L1 gives
||C_K||<=8sqrt(ell)=4/5, hence W_C<=16I/25. These are inequalities
for the full four-state loading map, not an unproved projection of it.
Its minimum-L2 source for target y is C_K^* W_C^-1 y and has squared
norm y^T W_C^-1 y. In particular

    (5/4)||y|| <= s_min(y) <= (2*10^6/3)||y||.

The sharper regular-light lower bounds are proved independently by the
exact symplectic dual kernels in notes/construction.md. The finite-band
optimization, its strict-budget smooth relaxation and finite-delta
error are in notes/optimization.md.

No least-cost source or band optimum is numerically computed here.
The validated computation supplies a coercivity bound; its role in the
analytic optimization is explicitly limited to that bound. These source
norms are linear response budgets, not nonlinear stress bounds, physical
energy estimates, a rolling mass gap, a cutoff or a UV-matching verdict.
