# S6.104: global relative-Proca first-order mean response

This child derives the finite state-dependent one-loop mean response of the
actual epsilon=1/200 constant-Proca clock. It does not solve the common
vacuum tadpole, an absolute semiclassical equation, or V/G/B UV matching.
Original P8 remains OPEN.

Retain S6.103's global seven-mode reference and its original physical
background a=(1+u^2)^2. Change only the Proca covariance at u=0 by

    Delta C(0,k)=eta b(k) I6,  eta=hbar lambda/kappa>=0,

where b is the fixed smooth nonnegative rotationally invariant band
supported in 1<|k|^2<4 and normalized in d^3k/(2pi)^3.
The added covariance is transported with the actual global Proca
generator on both legs. It is positive, has zero antisymmetric part,
and is smooth in spacetime. The target remains Hadamard. The scalar
and tensor fluctuation factors and all counterterm choices are unchanged.

The resulting physical vector density r and isotropic pressure s are
finite, obey -r/3<=s<=r, and satisfy r'+3H(r+s)=0. In the spatial
hat chart the lapse N is UNCHANGED, while a_physical=e(N,u) hat_a,
e=[(h-1+N^-2)/h]^(-1/4), h=(1+u^2)^3. The actual vector lapse
source is therefore F=r-3s/(2h), not r and not a four-dimensionally
conformally rescaled lapse source.

Let xi=delta log(hat_a), dp be the variation of
p=2 pi_trace/(3 hat_volume), and n=delta N. Fix the free matter
charge and zero initial xi,dp. With ell=1/[10(1+u^2)^6],
alpha=Hcal_pN, beta=ell(1-3/(2h)), and J=Jnew, the complete
linearized equations from the actual action are

    n=(alpha dp-3ell beta xi+F)/(2J),
    xi'=-dp/2+alpha n/3,
    dp'=-3H dp+ell beta n-3ell^2 xi+s.

The induced scalar field satisfies delta psi'=beta n-3ell xi,
with zero initial value. Its physical density and pressure both
change by -3ell^2 B, where B=xi+n/(2h) is the linear physical
log-scale response. Including the connection variation gives
the full homogeneous matter-plus-vector Ward identity at this order.

All six vector entries of S6.102's induced scalar density Hessian
agree independently with this mean constraint at u=0. The total
scalar-plus-vector density response there is at least 134r/135.
The earlier negative scalar SECOND VARIATION is not a negative
combined state response or negative full energy.

On |u|<=1/100, explicit bounds control every band momentum, the
two mean variables, lapse, physical scale, matter field and density.
An independent center stress-jet calculation bounds the proper
Hubble-derivative response by 5*10^11 eta.

A stronger all-time argument uses J>=1/[5(1+u^2)], positivity and
time-reflection symmetry of the ADDED covariance, and the weighted
phase (xi,100(1+u^2)^3 dp). It gives, for all real u,

    |xi| <=2.4*10^9 eta,
    |n| <=5.5*10^8 eta,
    |B| <=2.7*10^9 eta,
    |delta psi| <=4*10^8 eta.

The lapse tends to zero at both ends and xi has finite tail limits.
For 0<=eta<=10^-20 an explicitly defined exact spatial-frame
representative of this FIRST-ORDER profile is a smooth complete
bouncing metric. It is not an exact solution of the full field
equations; nonlinear and higher-loop remainders have not been
bounded by this statement.

The calculation is the difference of leading mean equations around
the same classical solution. State-independent subtractions and the
unchanged sectors cancel in this difference. No unknown absolute
reference source is set to zero, and no frozen tadpole is reassigned.
