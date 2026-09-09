# Differentiated subtraction and the new integrable tails

Use the same eighth-order comparison reference W=omega*S as
S6.59, not a new state: S=1+sum_(n=1)^4 P_n*t^n,
t=omega^-2. Its exact coefficient functions, bounds S>1/2 and
reference products are unchanged. Write R=D S and c=c1+R/(2S).
For the new exact-ODE-projected row the reference readout is
omega*F_j/(4a^3), where

    F_j=[t(rX-c*rY+c^2*rZ)+S^2*rZ]/S.

This is generally not the j-th time derivative of the approximate
mode readout. Projection by the exact ODE is essential.

The actual ordinary Proca AD4 bracket has coefficients
2, c1^2, P1^2+c1*B1-c1^2*P1 at orders t^0,t^1,t^2.
Here P1 is the original coefficient called P2 and
B1=D0 P1-2lambda P1. Apply
D0+[(1-2n)lambda-3H] successively to coefficient n to obtain
the actual differentiated subtraction AD_j. Its new energy
conservation identity with the old pressure subtraction is checked
at all three orders.

The unchanged general algebra identity is

    4t^2*S^3*(F_j-AD_j)
      =t^2{4t(rX-c1*rY+c1^2*rZ)S^2
            +2t(-rY+2c1*rZ)R*S+t*rZ*R^2
            +4rZ*S^4-4AD_j*S^3}.

It is used with new rows and new subtraction coefficients; no
parent readout is changed. For both polarizations and j<=5,
all coefficients t^0 through t^4 vanish exactly. All remaining
coefficients are bounded by their continuous nonnegative
polynomial envelopes, using |c1|<=2, t<=10^-6 and
1/(4S^3)<=2. Consequently

    |F_j-AD_j|<=T_j*t^3,
    |Q_projected,j-AD_density^(j)|<=T_j/(4a^3*omega^5).

The integer T_j for j=0,...,5 are:

    transverse:
    40976759,330164332,5423162383,407150071277,
    20454860860415,699227829249482

    longitudinal:
    39243120,343842778,6026341152,423718215792,
    25184980639455,1131096194226166.

A lower-reference failure control repeats the actual ordinary
energy calculation with W4. For the fourth time derivative of
longitudinal energy its numerator coefficient t^4 is 23808 at
u=0,z=1, not zero. This disproves that lower reference's separate
absolute-tail argument. It does not prove divergence of an exact
state derivative, and the old nonminimal j=2 failure control is
not silently reused.
