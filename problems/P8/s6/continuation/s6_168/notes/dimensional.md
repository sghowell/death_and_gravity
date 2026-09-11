# Restore the full dimensional finite parts

The two subtracted terms are a way to isolate a convergent
integral, not a choice that sets their finite parts to zero.
Use four-component Dirac trace and two occupied helicities,
d=3-2epsilon spatial dimensions, the full angular measure,
and MSbar at mu=mF. The same trace convention underlies
the frozen vacuum potential and two-point kernel.

The dimensional radial integral identity is

    integral d^d p/(2pi)^d (p^2)^a/(p^2+M^2)^b
    =(M^2)^(d/2+a-b)/(4pi)^(d/2)
      Gamma(a+d/2)Gamma(b-a-d/2)/[Gamma(d/2)Gamma(b)].

Meromorphic continuation gives

    -2 mu^(2epsilon) integral omega
      =M^4/Q [Ibar+3/2-log(M^2/mu^2)]+O(epsilon),

    Mdot^2/4 mu^(2epsilon) integral p^2/omega^5
      =Mdot^2/Q [Ibar-2/3-log(M^2/mu^2)]+O(epsilon),

with Q=16pi^2. After factoring out the pole, the two
MSbar-normalized analytic gamma factors are

    2 exp[(EulerGamma-ell)epsilon]Gamma(1+epsilon)
      /[(epsilon-1)(epsilon-2)],

    (1-2epsilon/3)exp[(EulerGamma-ell)epsilon]Gamma(1+epsilon),

where ell=log(M^2/mu^2). Their constant terms are one;
their first coefficients give 3/2-ell and -2/3-ell.
Dropping the d-dimensional angular factor would change
the derivative finite part. A cutoff finite term is
not substituted for this fixed prescription.

The restored finite energies per Dirac color/flavor are

    V_MS(M)=-M^4/Q [log(M^2/mF^2)-3/2],
    K_MS(M,Mdot)=-Mdot^2/Q [log(M^2/mF^2)+2/3].

At M=mF the derivative expression is exactly consistent
with the frozen two-point slope fF'(0)=4NY/(3Q):
the induced kinetic energy is -fF'(0) Phidot^2/2.

The relevant pole action is Ibar/Q times the integral of
sqrt(-g)[M^4-(partial M)^2], with the usual required
curvature counterterms. The external M is held as a
prescribed scalar under metric variation. In flat
homogeneous geometry, the 00 variation of R F(M) is
proportional to (g00 Box-partial0^2)F(M)=0. One cannot
drop those terms for pressures or curved geometry.

For comparison, [del Rio et al., v2](https://arxiv.org/pdf/1703.00908),
section IV equations44-48 with a=1, has the derivative
term p^2 sdot^2/[4(p^2+m^2)^(5/2)] and the fourth-order
s expansion of -2sqrt(p^2+(m+s)^2). Their section V
identifies the covariant counterterm structure. Here
the full M is retained and its dimensional finite parts
are derived above; no derivative expansion replaces
the exact-state remainder.

After the two subtractions, the three remainder
integrands decay at least as p^-2, p^-3 and p^-5.
They remain integrable for epsilon in a neighborhood
of zero. Hence the finite remainder can be evaluated
in three spatial dimensions before restoring these
explicit pole-subtracted local terms.
