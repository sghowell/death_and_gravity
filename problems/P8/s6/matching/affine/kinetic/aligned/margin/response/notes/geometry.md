# Exact point-chart lift of the first-order response

Do not identify the reduced v with the physical scale-factor
perturbation. Define a specific positive lift, on I only, by

    N_app=1+n,
    a_hat,app=a*exp(v),
    a_phys,app=a*exp[v+omega(u,N_app)],
    omega=-log[(h-1+N_app^-2)/h]/4.

At first order omega=delta*n, so this lift agrees with the
literal source calculation. The nonlinear functions above
are a chosen completion of the linear fields, not a solution
of the nonlinear equations.

For eta0<=10^-6, |n|<1/100. Thus N_app>0 and
X_app=N_app^-2 stays in [9/10,11/10]. The physical point map
is regular throughout the original timelike tube.

The exact derivatives are

    omega_N=1/[2N(1+(h-1)N²)],
    omega_u=h'(N^-2-1)/[4h(h-1+N^-2)].

For 99/100<=N<=101/100, h>=1 and |h'/h|<=3,
|omega_N|<1 and |omega_u|<3|n|. Here omega_u is the
partial derivative at fixed N and vanishes at N=1.
It follows that the physical log-scale shift z satisfies

    |z|=|v+omega(u,1+n)| <=225eta0.

This is below 1/2 on the entire stated source domain.
The elementary exp bound gives |exp(z)-1|<=2|z|, hence

    |a_phys,app/a-1| <=450eta0.

The actual Hubble rate of this lifted metric, in physical units,
is H_app=(H+z')/[tau(1+n)], not merely (H+v')/tau.
Using the previously proved bounds gives

    |z'-Hn| <49000eta0+4eta1,
    tau*|H_app-H/tau| <98000eta0+8eta1.

For eta0,eta1<=10^-6 this is less than 1/10.
The old dimensionless Hubble values at u=+-1/4 are
+-16/17, whose magnitudes exceed 9/10. The lifted Hubble
rate is therefore negative at the left point and positive
at the right. Positivity of N_app and continuity imply an
interior local minimum of the lifted physical scale factor.
Neither uniqueness nor a quantitative positive second
derivative at that minimum is asserted.

For the actual vector example both source norms are below
10^-14. The rounded bounds are

    |N_app-1| <9.4*10^-13,
    |a_phys,app/a-1| <4.5*10^-12,
    tau*|H_app-H/tau| <10^-9.

These are bounds on the specified fixed-source response and
its chart lift. They are not a uniform error bound to an
unknown self-consistent quantum solution, do not establish
global quantum completeness outside I, and do not bound the
quantum-corrected perturbation symbol.
