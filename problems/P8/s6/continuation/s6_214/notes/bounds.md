# Full all-momentum form norm and exact finite-band response

For a possibly complex Fourier direction D=(n,beta,Q), with Q complex symmetric, define

    sigma_D=|n|+(5/2)||Q||op+a||beta||2.

The trace inequality |trQ|<=3||Q||op and the triangle inequality give ||B_Q||op<=5||Q||op/2. This is a safe bound, not asserted optimal. The scalar constraint first form is bounded by3||Q||op/2. Together with the sharp shift bound and lapse identity form,

    ||S_D||op<=sigma_D.

Submultiplicativity controls the full symmetrized spatial product by25||Q_D||op||Q_G||op/4. The constraint second product is smaller. Lapse-spatial terms add

    (5/2)(|n_D||||Q_G||op+|n_G||||Q_D||op).

The full sum is bounded by sigma_D sigma_G. The difference is an explicit polynomial with nonnegative coefficients in all direction norm variables. No reality of a single Fourier mode is used.

Because F_k^tF_k=M0(k), the matrix F_k M0(k)^(-1/2) is an isometry on the full six-dimensional canonical space. Hence

    ||M0(k)^(-1/2) M_D(k,q) M0(q)^(-1/2)||op<=sigma_D,

and likewise for the second form with sigma_D sigma_G, for all k,q including zero and arbitrarily separated momenta. This is not a same-momentum estimate. The same actual symplectic normalizers convert it to sqrt(omega_k omega_q) times the direction form bound.

The exact S195 off-diagonal covariance tangent and its full current/contact formula apply unchanged to these full ADM matrices. The actual prepared covariance, both source terms and both distinct propagators remain. Its balanced covariance and propagator-product bounds are4e12 and4e9, so the same source1e23, readout1e24 and contact3e13 constants hold with sigma replacing the tracefree direction norm.

At the same sharp computational band |k|,|q|<=K, K>=1000,

    memory <1e24 K^5 sigma_D(t) integral_left^t sigma_G(s)ds,
    contact <1e13 K^4 sigma_D(t)sigma_G(t).

The contact has only the one-leg internal mask; empty propagation overlap at large external P does not justify deleting nonzero spatial/lapse contacts. At K=1e16 on the unit slab, the total is below2e104 for unit directions. These are exact finite-regulator statements with growing powers, not a full ADM ultraviolet subtraction or a physical cutoff.

The unreduced sigma norm is explicitly a metric-coordinate norm. No reduced canonical scalar/metric normalization, inverse, stability or quantum-background estimate follows from it.
