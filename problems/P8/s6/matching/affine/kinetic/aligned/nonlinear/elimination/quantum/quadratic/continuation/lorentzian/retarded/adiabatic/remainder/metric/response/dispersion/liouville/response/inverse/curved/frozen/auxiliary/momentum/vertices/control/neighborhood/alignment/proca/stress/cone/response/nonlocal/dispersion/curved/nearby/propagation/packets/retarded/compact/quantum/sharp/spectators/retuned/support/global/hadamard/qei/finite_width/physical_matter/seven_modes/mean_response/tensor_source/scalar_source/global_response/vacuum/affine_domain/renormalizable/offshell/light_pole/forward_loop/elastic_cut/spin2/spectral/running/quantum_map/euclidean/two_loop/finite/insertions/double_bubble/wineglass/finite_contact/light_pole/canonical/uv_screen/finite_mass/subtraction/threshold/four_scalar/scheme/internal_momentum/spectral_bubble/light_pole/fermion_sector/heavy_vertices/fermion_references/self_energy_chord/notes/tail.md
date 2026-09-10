# Soft projection, exact radial integral and coefficient bound

Work first at a common regulator, with every proper
self-energy and its MS mass/kinetic counterterms
paired. The complete 24-word subset in each sector
is Lorentz invariant and closed under S4 external
permutations. Its soft degree zero is a quartic
contact. Odd degrees vanish. Any degree-two
Lorentz scalar is a linear combination of external
Gram products. S4 averaging leaves only the sums
of diagonal and off-diagonal products; momentum
conservation relates them. It is therefore
proportional to sum p_i^2, constant on the
mass-one shell. None of these degrees contributes
to the forward quadratic coefficient b2.

The remaining overall divergence is precisely the
degree-zero local quartic term. Its regulated
subtraction is made before claiming a continuum
amplitude. Changing that contact cannot tune b2.
Higher-degree coefficients of the properly paired
subset are ultraviolet integrable; the following
pointwise bounds justify interchanging their
projection and the outer integration. The symmetry
argument is for the complete regulated subset,
not an assumed permutation symmetry of an
individual unshifted loop integrand.

Each marked box has five fermion propagators.
The four-dimensional Dirac trace bound gives

    |integrand| <=128 N Y^2 B
                   [log(1+q^2/mF^2)+6]/S^2,

where B=(3Y+12a Cf)/Q and N=6.
For the pointwise analytic soft-scaling series,
Cauchy's bound on radius R(q)=sqrt(S)/360 gives
the degree-four-and-higher tail at zeta=1 as at most

    2 times the integrand majorant times R(q)^(-4),

because R(q)>=2. This projection is performed
BEFORE integrating a potentially divergent
unprojected expression. It leaves the integrable
majorant

    256*360^4 N Y^2 B
       [log(1+q^2/mF^2)+6]/S^4.

The four-dimensional radial measure is
(1/Q) y dy, y=q^2. Scaling y=mF^2 t gives

    integral_0^infinity y[log(1+y/mF^2)+6]/
                        (mF^2+y)^4 dy
      =41/(36 mF^4).

Indeed u=1/(1+t) turns the two moments into
integral_0^1 u(1-u) du=1/6 and
integral_0^1 -u(1-u)log(u) du=5/36.
Both are checked exactly and independently
by numerical quadrature.

Include the 24 words. The exact prefactor is

    24*256*360^4*(41/36)=117528330240000<2*10^14.

The uniform projected-amplitude bound is therefore

    E <=2*10^14 N Y^2(3Y+16a)/(Q^2 mF^4).

All estimates hold uniformly on the unit forward
disc with an open-neighborhood margin. Cauchy's
formula in s then gives the same E for |b2|.
This is a second Cauchy step, after the convergent
outer integral has been established.

The actual rational upper bound is approximately
6.3285927498508 times 10^-1404, relative to the
tree 4 lambda approximately
1.5821481874627 times 10^-804. The API uses
mF>=720 and a conservative 0<Qlo<=144<16 pi^2.
These bounds cover both scalar and gauge
self-energy chords, not the other 36 words
in either primitive sector. A positive error
majorant is not a sign claim for the amplitude.
