# The isotropic dimensional factor changes the finite answer

Keep d=3-2epsilon before subtraction. The isotropic
pressure operator uses1/d, not a fixed1/3. With the
same four-component trace and MSbar convention,

    J0=integral p^2/omega^5
       =4/Q[Ibar-log(M^2/mu^2)-2/3]+O(epsilon),
    J1=integral p^4/omega^7
       =4/Q[Ibar-log(M^2/mu^2)-16/15]+O(epsilon),
    J2=integral p^2/omega^7=8/(5Q M^2),

where integrals have the full d-dimensional measure.
J2 is finite at epsilon=0. The full pressure term is

    -M Mddot J0/(2d)+Mdot^2 J1/(4d)
        +3M^2 Mdot^2 J2/(2d).

The normalized gamma factors for J0,J1 are

    (1-2epsilon/3)C(epsilon),
    (1-2epsilon/3)(1-2epsilon/5)C(epsilon),

    C=exp[(EulerGamma-ell)epsilon]Gamma(1+epsilon),
    ell=log(M^2/mu^2).

The epsilon term of1/d cancels the -2/3 constant in
the acceleration coefficient and changes the other
finite coefficient. The restored MS pressure is

    P0_MS=-V_MS(M),
    P2_MS=[2M Mddot ell/3-Mdot^2 ell/3
              +2Mdot^2/3]/Q.

Dimensional integration by parts establishes
P0=-rho0 before finite subtraction. For comparison,
[del Rio et al.](https://arxiv.org/pdf/1703.00908),
v2, equations56-57 at a=1 reproduce the Mddot and
Mdot^2 coefficients when the full-M expressions
are expanded to the indicated source orders.

The inherited energy term is
rho2_MS=-Mdot^2(ell+2/3)/Q. Their difference is

    P2_MS-rho2_MS=
      partial_t^2[M^2(log(M^2/mu^2)-1)]/(3Q).

It is the flat local tensor improvement

    -(eta_mu_nu Box-partial_mu partial_nu)F/(3Q),
    F=M^2(log(M^2/mu^2)-1).

Its 00 component vanishes for homogeneous F; each
diagonal spatial component equals Fddot/(3Q).
Zero background curvature is not a reason to drop
this metric variation. The tensor formula fixes
the convention without guessing a curvature-action
sign from a different signature convention.

The required curvature poles remain in the common
MS prescription; additional finite curvature terms
are set to zero here. This is an explicit extension
of the previous energy-only prescription, not a
claim that flat T00 fixed those finite terms.
Terms invisible to this flat tensor do not establish
the physical curved Newton-reference dictionary.
