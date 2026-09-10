# Project first, then integrate both loop momenta

At fixed real q,l the analytic word is bounded
on the joint zeta circle by

 M(q,l)=256 N Y^2 (Y+4aCf) /
         [S_q^(3/2) S_l^(3/2) (q-l)^2].

Cauchy's Taylor estimate and R>=2 imply

 |sum_{n>=4} a_n| <= M/[R^4(1-1/R)]
                  <=2 M/R^4.

Perform this soft projection before either
UV integral. After angular integration the
positive double radial majorant, apart from
512 times 360^4 N Y^2 (Y+4aCf)/Q^2, is

 integral_0^infinity dx dy x y /
 [max(x,y)(m^2+x)^(3/2)(m^2+y)^(3/2)
  min(m^2+x,m^2+y)^2].

Here Q=16 pi^2, since each normalized
four-dimensional loop measure is x dx/Q.
In x>=y the x integral is
integral_y^infinity dx/(m^2+x)^(3/2)
=2/sqrt(m^2+y). The two equal radial regions
therefore give

 4 integral_0^infinity y dy/(m^2+y)^4
   =2/(3m^4).

For the equality use y=m^2 t and u=1/(1+t);
the remaining compact moment is
integral_0^1 u(1-u)du=1/6. This proves both
UV and massless-chord IR integrability, not
merely numerical smallness below a cutoff.

For any positive soft degree n the same
argument gives 16/[n(n+2)m^n]. Thus every
positive-degree coefficient is absolutely
integrable, including degree two. The degree
zero contact is subtracted in a common
translation/rotation-invariant regulator
before removing it. The complete twelve-word
subset is S4 invariant. Its integrated odd
soft degrees vanish by Lorentz invariance.
A degree-two S4 Gram polynomial has only
the diagonal and off-diagonal orbits; momentum
conservation relates them. It is proportional
to sum p_i^2 and is constant on the mass-one
shell. Consequently degrees below four have
zero forward b2. This argument concerns the
complete regulated subset, not an individual
unshifted integrand.

The integrated projected amplitude is uniformly
holomorphic in the unit forward disc. A second
Cauchy coefficient estimate gives the same
bound on b2. Summing twelve words yields

 E_opposite <= 68797071360000
      N Y^2 (Y+4aCf)/(Q^2 m^4)
    < 10^14 N Y^2 (Y+16a/3)/(Q^2 m^4).

Use the shared rational upper bounds on Y,a,
Q>=144, and m=10^200. Scalar/gauge pieces
are kept separately. The combined rounded
bound is approximately 1.05476545830846e-1404,
relative to 4 lambda approximately
2.63691364577116e-805. These positive majorants
do not establish an amplitude sign.

A first-order finite field or parameter
conversion on this new order-two family is
order three. Order-two conversions multiplying
old order-one amplitudes remain separate
ledger entries. No such entry is discarded.
