# Full metric vertices and their energy norm

Let G,D be arbitrary real symmetric tracefree shear
directions and C=(GD+DG)/2. Differentiating the exact
Hamiltonian at gamma=0 gives the symmetric matrices

    M_G=diag(-a m^2 G+a^-1 C_k^T G C_k, a^-1 G),
    M_GD=diag(a m^2 C+a^-1 C_k^T C C_k, a^-1 C).

The kk^T temporal-constraint term is shear independent
in this unimodular chart, not absent from the dynamics.
Independent differentiation of both inverse spatial
metric factors in the magnetic action gives the same
first and second vertices, including the mixed terms.

At the base, physical frame fields are E_phys=pi/a^2,
B_phys=(k cross A)/a^2 and massA=mA/a. For a tracefree
direction, the pressure trace and temporal mass trace
drop from T_G, leaving
T_G=-E_phys^T G E_phys-B_phys^T G B_phys+massA^T G massA.
Consequently H_G=-a^3 T_G/2. The external coordinate
current is J_G=-H_G. Away from the base its readout
uses the full derivative of exp gamma, not the base
physical component silently kept fixed.

The second vertex is generally nonzero. In particular
M_GG is positive semidefinite and is nonzero when G
is nonzero. In a positive full-rank finite-mode covariance
its expectation is strictly positive. Removing its
response contact therefore changes the answer.

For every n>=1, the nth derivative of exp gamma at0
is the average of all ordered products of its n
directions. The norm is at most their norm product.
It is self-adjoint for symmetric directions. The
positive mass/magnetic/electric decomposition then gives

    ||M0^-1/2 D^n M_0[G1,...,Gn] M0^-1/2||op
       <= product_j ||Gj||op.

The constant constrained term only increases the
reference denominator. Noncommuting directions are
included. For canonical h=sqrt(kappa)gamma/2 the
nth vertex has the additional factor(2/sqrt(kappa))^n.
This is a tensor/shear chart statement, not a reduced
scalar-metric cutoff estimate.
