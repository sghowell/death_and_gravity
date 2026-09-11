# Complete mixing-series tail and a uniform transition-integral bound

Write I1(p)=integral v(t) exp(-i theta(t)) dt for
the first transition-coupling Dyson term. The phase
contains the full varying omega; it is not a
mass-amplitude or Feynman-loop expansion.

Since |M'|=Delta s'(t/tau)/tau and omega>=
E=sqrt(p^2+m0^2),

    eta=integral |v| dt <=p Delta/E^2
        <=Delta/(2m0)<1/100.

Every interaction generator is off diagonal.
Only odd ordered products contribute to beta.
Absolute convergence of the Dyson series gives

    |beta-I1|<=sinh(eta)-eta
              <=eta^3 exp(eta)/6
              <=eta^3/[6(1-eta)]
              <=eta^3/5.

The ordered integration volumes supply the
factorials. The last coefficient follows from
eta<=1/100 and 50/297<1/5. This retains ALL
higher transition terms of the quadratic
evolution, not all quantum Feynman loops.

Two integrations by parts have zero boundary
terms. Up to the irrelevant overall sign,
their integrand is

    L2v=partial_t[(partial_t(v/(2omega)))/(2omega)]
       =p/8[M'''/omega^4
             -7M''omega'/omega^5
             -3M'omega''/omega^5
             +15M'(omega')^2/omega^6].

For the full frequency,

    |omega'|<=|M'|,
    |omega''|<=|M''|+|M'|^2/omega.

Therefore

    |I1|<=p/8 integral[
        |M'''|/E^4
        +10|M'M''|/E^5
        +18|M'|^3/E^6]dt.

The three profile integrals give, with r0=Delta/m0,

    |I1|<=p Delta/(8tau^2 E^4)
             [36+10r0+36r0^2]
          <5p Delta/(tau^2 E^4).

Combining this with the complete higher-order
transition tail gives the uniform bound

    |beta_p|<=A(p)+B(p),
    A=5p Delta/(tau^2 E^4),
    B=p^3 Delta^3/(5E^6).

Both terms decay as p^-3. Mode unitarity
separately gives |beta_p|<=1; no sign is
assumed for the transition integral.
