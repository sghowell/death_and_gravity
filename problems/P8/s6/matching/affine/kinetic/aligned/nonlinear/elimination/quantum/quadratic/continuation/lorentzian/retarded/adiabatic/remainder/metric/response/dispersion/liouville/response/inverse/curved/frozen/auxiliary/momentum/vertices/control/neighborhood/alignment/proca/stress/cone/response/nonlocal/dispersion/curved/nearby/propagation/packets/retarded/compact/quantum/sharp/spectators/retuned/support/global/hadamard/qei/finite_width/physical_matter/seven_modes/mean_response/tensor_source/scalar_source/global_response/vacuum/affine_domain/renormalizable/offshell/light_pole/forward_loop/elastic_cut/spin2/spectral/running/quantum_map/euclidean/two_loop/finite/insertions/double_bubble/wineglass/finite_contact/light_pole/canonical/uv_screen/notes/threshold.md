# Leading heavy-fermion gauge threshold

Take constant Phi with |y Phi|<mF, so both flavor masses
mF+y Phi and mF-y Phi are positive. Their squared-mass
determinant is (mF^2-y^2 Phi^2)^2. Subtracting its vacuum
logarithm gives 2 log(1-z), z=y^2 Phi^2/mF^2.

For a Dirac fundamental the leading heavy threshold is
a F^a F^a log(det Mdagger M)/(96 pi^2), using the
conventions and sourced low-energy theorem in the literature
note. An independent normalization check starts from the
one-Dirac contribution delta b=2/3 to gauge running.
In -Fcal^2/(4g_s^2), threshold running supplies the coefficient
delta b log(m)/(32 pi^2). Restoring Fcal=g_s F and adding
the two masses reproduces the determinant expression.

The first Phi derivative vanishes. The coefficient of
Phi^2 F^a F^a is instead

    C = -a y^2/(48 pi^2 mF^2),

strictly negative when the two couplings are nonzero.
Opposite Yukawa signs do not cancel an even threshold.

The code differentiates the literal determinant and the
reference expression independently; it also checks the
first eight field-series coefficients. For 0<=z<1,

    |2 log(1-z)+2z|
      = 2 sum_{n>=2} z^n/n <= z^2/(1-z).

The nonnegative difference between the bound and the exact
magnitude vanishes at zero and has derivative z^2/(1-z)^2.
At z<=1/4 the remainder is at most 4z^2/3.

This is a constant-field and leading derivative expansion.
It is not a bound on the finite-mass momentum form factor.
In particular, scalar external on-shell p_i^2=1 must not be
replaced by zero without a derivative remainder calculation.
The next note is an exact calculation in the leading
matched local operator theory, not that missing remainder.
