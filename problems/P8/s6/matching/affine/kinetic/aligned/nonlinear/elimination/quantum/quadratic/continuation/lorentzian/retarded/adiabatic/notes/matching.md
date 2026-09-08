# Radial finite term and the dimensional counterterm

For a polynomial c(z,D), the physical radial integral of
omega^(1-2j)c/4 has the coefficient

    G_j(D)=2sqrt(pi) sum_l c_l(D)(D/2)_l/Gamma(j+l-1/2),

relative to m^(4-2j)/(64pi^2), with Gamma(j-2+epsilon) restored
and D=3-2epsilon. Therefore

    P_j=(-1)^(2-j)G_j(3)/(2-j)!,
    E_j=-2(-1)^(2-j)G_j'(3)/(2-j)!,
    finite_radial(mu=m)=Harmonic(2-j)P_j+E_j.

The frozen S6.52 radial formula is reused. The varied mass-readout
polynomials, including their D dependence, are new.

On the actual clock with alpha=4/[9(1+u^2)^3] and
beta=28/[81(1+u^2)^3], all three poles satisfy exactly

    m^(4-2j)P_j=2 E_loop,2j

for the normalized homogeneous Euler operator of S6.63.
The factor two converts its 1/(32pi^2) convention to the current
1/(64pi^2) convention. This is checked through the full curved
two- and four-derivative terms, not only in flat space.

The bare-action counterterm contributes +2 partial_D E_loop,D
relative to 1/(32pi^2). Thus its contribution in the present
normalization is +4 partial_D E_loop,D. It is added to the radial
finite part only after the common-dimensional calculation.
The existing zero-derivative potential continuation is D-independent
and adds no such normalized mass-source evanescent operator.

The result is the following compact finite local action component:

    (64pi^2)^(-1) integral du a^3 [A n''^2+B n'^2+C n^2].

Put r=1+u^2. Separating adiabatic orders 0,2,4 gives

    A0=B0=0, C0=-1652 m^4/(6561 r^6),

    A2=0, B2=-800 m^2/(6561 r^6),
    C2=-128 m^2(917u^2+164)/(6561 r^8),

    A4=628/(6561 r^6),
    B4=80(1483u^2+459)/(19683 r^8),
    C4=64(3651u^4+4285u^2+338)/(6561 r^10).

Varying these compact actions exactly reproduces the directly
computed finite response. This self-adjointness is a nontrivial
check: the radial finite part without the evanescent counterterm
fails it. Its next-to-highest derivative coefficient has a
nonzero adjoint defect. The zero-derivative finite response
also equals minus the second mass-source derivative of the
frozen full finite potential, so it must not be added a second time.

The displayed term is the local subtraction component, not a
complete finite physical response by itself. Additional finite
Y^2 Wilson coefficients are not fixed here. The nonlocal
subtracted exact response, its common regulator limit and its
quantitative bounds remain open.
