# Actual uniform covariance remainder and infinite tail

Use the literal background a_max=25/16 and define

    nu_minus^2=m^2+(99/100)|k|^2/a_max^2.

Since ||gamma||<=1/100,
exp(-gamma)>=exp(-1/100)I>=(99/100)I. Hence
omega>=nu_minus>=m. Also exp(-gamma)<2I and a>=1
give omega<=3nu_minus. The initial frequency nu is
at least nu_minus. These bounds are uniform in k,
the direction of k and every admitted shear history.

For nu_minus>=K=1e16, combine state.md and reference.md:

    ||r_actual-rhat||
      <1e6(1e33/K+C) nu_minus^-10
      <1e30 nu_minus^-10.

This is an actual-state error over the complete unit
slab, not a finite set of sampled momenta and not
a change to the Cauchy data. The actual graph and
the reference remain inside the radius1/10 covariance
ball, because ||rhat||<1/100 and the displayed error
is tiny in this high band. Therefore frame.md gives

    ||Sigma_actual-Sigma(rhat)||op
      <2e31 nu_minus^-10.

Sigma is the full real symmetrized covariance of
(Q_1,Q_2,Q_3,P_1,P_2,P_3), including correlations.
The bound applies to finite admitted prescribed shear
histories, without assuming commuting time directions.
It is not yet an estimate for derivatives with respect
to the shear amplitude parameter.

The standard three-dimensional continuum measure is
d^3k/(2pi)^3. Put b=(99/100)/(25/16)^2. The radial
Jacobian satisfies

    k^2 dk =b^-3/2 nu_minus sqrt(nu_minus^2-m^2) dnu_minus
            <=b^-3/2 nu_minus^2 dnu_minus,
    b^-3/2<4.

Using omega<=3nu_minus and pi^2>9 proves

    integral_(nu_minus>=K) omega
       ||Sigma_actual-Sigma(rhat)||op d^3k/(2pi)^3
      <[3(2e31)4/(12*9)] K^-6
      <1e-65.

The infinite radial integral is evaluated, not replaced
by a finite numerical momentum grid. K partitions a proof;
it is not a physical cutoff and does not remove modes.
The finite low band still requires its own estimate.
The reference covariance itself retains UV-divergent
terms that must be subtracted and matched to the original
covariant prescription. The present absolutely integrable
difference does not make that matching automatic.

No small inverse, finite-coupling feedback remainder,
nonlinear quantum background or full parent matching
follows from this positive covariance-tail integral.
