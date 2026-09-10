# Finite Cauchy subtraction and local cut removal

Fix T=6, in invariant mass-squared units with the light
mass one. Initially off the real cut, define

    C_rho(s)=(1/pi) integral_0^6 rho(x)/(x-s) dx,
    F_low(s)=C_rho(s)+C_rho(4-s).

This is a specified light-cut subtraction, not a change
of the parent action or a finite counterterm fitted to
a desired sign. No infinite-energy representation is
assumed in defining these finite integrals.

On an interior cut point the elementary boundary
identity for 1/(x-s-i0) gives
C_rho(s+i0)-C_rho(s-i0)=2i rho(s).
The crossed u=4-s has the opposite boundary sign.
For 0<s<4 the combined jump is therefore

    2i[rho(s)-rho(4-s)].

It is exactly the first two-gauge jump computed in
S6.129, not its leading approximation. On the real
part of |s-2|<1 there are no endpoints of these
Cauchy integrals. Subtracting F_low from the two
boundary germs of that channel removes its known
jump. Their locally bounded continuations glue
holomorphically by the contour form of Morera's
theorem. This is a local first-channel, fixed-order
statement; it does not establish all-loop analyticity
of the candidate's full scattering amplitude.

## Removable representation and branch convention

Separate the numerator before differentiating:

    C_rho(s)=(1/pi) [
      integral_0^6 (rho(x)-rho(s))/(x-s) dx
      +rho(s)(Log(6-s)-Log(-s)) ].

The divided difference is removable at x=s.
On the unit disc about s=2, the local continuation
from the upper physical boundary uses
Log(-s)=Log(s)-i pi. The crossed term uses
Log(-u)=Log(u)+i pi. The lower germ interchanges
these two signs. These are analytic local germs,
not a declaration that the original physical
upper and lower boundaries were identical.

The difference is crossing odd about s=2.
All even center derivatives of the two germs
therefore agree, even for the finite-mass rho.
This does not recover the original gapped proof:
a neighborhood boundary mismatch was present
before the explicit subtraction.

For the leading rho=pi K x^2, the primitive is

    integral x^2/(x-s) dx
      =x^2/2+s x+s^2 Log(x-s).

At general T>4 its center second coefficient is

    b2_low^0/K =
      2 log((T-2)/2)-8/(T-2)-4/(T-2)^2-3.

At T=6 this is 2log(2)-21/4.
The derivative with respect to T is
2T^2/(T-2)^3, independently obtained by
differentiating the finite-integral endpoint kernel.
This records the subtraction-scale dependence;
the cutoff is not silently optimized or absorbed
into a renormalization of the physical action.
