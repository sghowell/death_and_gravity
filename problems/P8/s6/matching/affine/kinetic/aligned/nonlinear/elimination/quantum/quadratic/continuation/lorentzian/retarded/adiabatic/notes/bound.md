# Continuous bound for the finite local source operator

The compact local action has normalized Euler operator

    O n=2C n-2(D_u+3H)(B n')
        +2(D_u+3H)^2(A n'').

For these exact coefficients, a^3 A and a^3 B2 are constant.
Hence the third-derivative coefficient vanishes, as does the
first-derivative coefficient at adiabatic order two. The full
Euler operator, rather than just the compact action coefficients,
is used in the estimate.

On I=[-1/2,1/2], every normalized coefficient is a polynomial in u
divided by a positive constant times (1+u^2)^p. Bound the numerator
by the sum of absolute polynomial coefficients times (1/2)^degree
and the denominator below by its positive constant. This is a
continuous interval bound, not sampling at a finite grid.

After removing m^(4-2j), the upper bounds for the coefficients of
n,n',n'',n''',n'''' are

    order 0: (3304/6561,0,0,0,0),
    order 2: (100672/6561,0,1600/6561,0,0),
    order 4: (23288/729,149720/19683,132760/19683,0,1256/6561).

Let ||n||_C4 be the maximum over those five time derivatives on I.
Since pi^2>9, divide the mass-weighted sum of these bounds by
576(M tau)^2. At M tau=10^24 and m0 tau=1000 the exact result is

    206506392019063 /
    236196000000000000000000000000000000000000000000000000

and is strictly less than 10^-39. Thus the normalized finite
local component maps C4 source data into C0 outputs with norm
below that value on the stated interval.

This estimate does not cover the nonlocal subtracted response,
spatially varying sources or a coupled inverse/gravitational
feedback map. It cannot by itself establish a contraction,
quantum stability, propagation cones, interaction cutoff or
UV matching. Higher-derivative local terms are not resummed
into exact extra modes of the finite EFT.
