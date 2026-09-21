# Exact common-basis mass moments

Let M=1+(n-1)z. Add the S344 primitive curvature weight only at j0:

    H0=K0-4z^2(1-z)^5/15
      =-2(z-1)^4(151z^3+44z^2+12z+3)/105,
    Hj=Kj for j=1,2,3,
    Jj(n)=int_0^1 Hj(z)/M^(4-j) dz.

Set y=M. A weight of degree d gives a polynomial numerator divided by
(n-1)^(d+1)y^(4-j). Integrate every monomial exactly, using log(y)
for exponent-1. Verify the full primitive derivative and both endpoints.
The source computes all four rational-log closed forms, their removable
equal-mass limits and exact rational mass calibrations.

The box is the already matched S345 integral

    c_box=int_0^1 [-4t(18+18t+92t^2+363t^3)/315]/M^5 dz,
    t=z(1-z).

Its independently integrated primitive agrees with the complete frozen
box expression. This is the same comparison basis, not a new box term.
All moments concern positive internal masses and analytic-origin jets;
they do not justify a Taylor approximation on the physical timelike cut.
