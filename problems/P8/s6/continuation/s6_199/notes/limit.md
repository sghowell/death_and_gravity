# Controlled spectral tail, not a physical cutoff

On x=m^2/s in [0,1/4], beta<=1. The polynomial bounds follow from

    30-(13+56x+48x^2)=(1/4-x)(48x+68)>=0,
    1-(1-4x+12x^2)=4x(1-3x)>=0.

Consequently

    0<=rho2(s)<=s^2/(128 pi^2),
    0<=rho0(s)<=s^2/(384 pi^2).

For a computational squared-spectral-mass limit Lambda>=4m^2 and
|z|<=2m^2, the discarded part obeys

    |D2_tail|<=2|z|^3/(128 pi^2 Lambda),
    |D0_tail|<=2|z|^3/(384 pi^2 Lambda),

using |s-z|>=s/2 and integral_Lambda^infinity s^-2 ds=1/Lambda.
The bounds control the actual improper dispersive integrals, not a
sampled finite interval. The tests independently integrate the tail
from v_Lambda=sqrt(1-4m^2/Lambda) to 1 for three different limits.

Lambda is not the physical effective-theory cutoff and no high-energy
degrees of freedom are removed. Positivity of this one Gaussian
external-metric cut does not prove a dispersion relation, contour bound,
omitted-order estimate or optical theorem for the full interacting
parent. The all-energy spectral density remains unbounded as s^2
despite multiplication by 4/kappa. No order reduction, inverse pole
removal or global scattering-history assumption is introduced.
