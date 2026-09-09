# An actual all-order state with an explicit initial envelope

Use the full S6.99 canonical Riccati recursion, not a leading vacuum
replacement. Exact finite time jets over QQ(sqrt(17985),i) solve the
actual recurrence through order eight, checking every residual
through its declared triangular time order. Their even initial
coefficient norm bounds are

    ||R2(0)||<32, ||R4(0)||<506,
    ||R6(0)||<9000, ||R8(0)||<600000.

These bounds use a rational isolating interval for sqrt(17985);
the large cancelling algebraic terms are not evaluated as rounded
decimals. An independent degree-nine exact calculation is also
contained in all 196 relevant whole-interval ball-jet entries.

The complete scaled Hamiltonian has even diagonal h0,h2 and odd
mixed h1. The positive-root leading graph is even in u.
The frequency-sum recursion therefore proves R_n(-u)=(-1)^n R_n(u)
by induction. It also gives imaginary even coefficients and real
odd coefficients. In particular EVERY odd R_n(0) vanishes, not
merely the finitely tested ones.

Choose real Borel cutoffs for n<=8 equal to one at k>=16.
For n>8 choose their radii successively large enough for all of
S6.99's symbol seminorm requirements AND

    ||cutoff_n(k) k^(-n) R_n(0)|| <= 2^(-n) k^(-8).

This is possible because n>8 and the initial coefficient is finite.
These added requirements do not disturb the all-order Borel
construction on a compact neighborhood. Thus, with
R6_trunc=sum_(n=0)^6 k^(-n) R_n, the resulting genuine symbol obeys

    A_B(0,k)=0,
    ||R_B(0,k)-R6_trunc(0,k)|| <= 600001 k^(-8), k>=16.

The exact initial leading B0 has 0.9 I<=B0 and ||B0||<12.8.
The sum 32/16^2+506/16^4+9000/16^6+600001/16^8 is <0.15.
Hence the ALL-ORDER initial B_B=-Im R_B has

    (3/4)I < B_B < 13 I, k>=16.

Take its positive Gram covariance in the original normalized
canonical chart. Below k=16 use the fixed B0 graph in the regular
normalized chart. Smoothly mix the two covariances with a radial
cutoff equal to zero up to 16 and one from 32 onward. This convex
mixing preserves positivity and the same exact CCR; it need not
preserve purity. The zero-momentum end is smooth in the spatial
momentum vector. The change at finite momentum is smoothing.

At u=0 the normalized chart transition is diagonal with powers
r^(-3/2),r^(1/2),r^(3/2),r^(-1/2), where r=mu/k.
For k>=16 its squared norm is <1.01. The full initial regular
covariance, with hbar/(2 kappa) removed, has norm <16.
Also ||V(0)^-1 B_B(0,k)^(-1/2)||<2.

Transport this initial covariance with the exact original GLOBAL
density evolution. S6.99's all-order signed-wavefront argument gives
a genuine global generalized Hadamard state; the finite-momentum
modification is smooth. No finite truncation is called Hadamard.
The bounds hold for the family of Borel choices just specified;
no arbitrary numerical reference functional is assigned.
