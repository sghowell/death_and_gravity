# A finite interpolation polynomial with two exact flat jets

Let T8 be the polynomial integral

    T8(X)=51480 integral_0^X t^7(1-t)^7 dt
         =6435X^8-40040X^9+108108X^10-163800X^11
          +150150X^12-83160X^13+25740X^14-3432X^15,
    q8(X)=1-T8(X).

The normalization gives T8(1)=1 and
T8(X)+T8(1-X)=1. Thus q8-1 is divisible by X^8,
while q8 is divisible by (1-X)^8. At the vacuum,
q8(0)=1 and derivatives one through seven vanish.
At the clock, q8 and all derivatives through seven
vanish. The leading coefficients are

    q8(X)-1=-6435 X^8+O(X^9),
    q8(X)=6435 (X-1)^8+O((X-1)^9).

This is a finite polynomial, not the separate rational
step in the actual analytic target. It introduces no
denominator or infinite field series into the map.
No uniform bound on q8 outside the two stated domains
or global inverse is claimed.

Use X=(partial Psi)^2/kappa in fixed light-mass units.
Keeping reference mass factors explicit, the dimensionless
D-dimensional argument is

    X_D=mu^(2epsilon)(partial Psi_D)^2/(kappa m_Phi^4),

with m_Phi=1 in the numerical formulas and dimensionless
kappa. Each X_D^j supplies mu^(2j epsilon). Use the
same cubic R_D, with its restored light-mass factors,
as in S6.160, and apply this common lift to the action,
counterterms, Jacobian and physical source before
finite-part operations.
