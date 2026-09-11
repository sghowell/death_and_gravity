# Full causal fixed-metric light-force bound

The conditional action in the independent canonical variable A is

    -F(A)^2/4+A^2/(2zeta)-A.J+kappa S^2/2,
    J=sqrt(kappa/zeta)S.

At fixed metric, varying u in a compact direction eta BEFORE
eliminating the vector gives

    delta S_vector = kappa integral sqrt(-g) (S-W).D S[eta].

Both the temporal constraint and source contact are present.
Take its expectation in the specified sourced state: the centered
Gaussian mean vanishes, so W is replaced by the actual RETARDED
mean Wbar. The fixed-metric connected determinant has zero
u variation in this specified conditional measure. This is an
in-in force expectation obtained from the unintegrated action,
NOT the derivative of one half S Gret S. The latter would
symmetrize Gret with Gadv and introduce a future dependence.

Classically, the scalar variation at fixed original affine trace
differs from the independent W-coordinate variation by the vector
Euler equation times the field-map variation; that equation
vanishes at the retarded mean. No full affine quantum Jacobian
or measure is thereby computed.

The mass-energy term controls the spatial mean, and the ACTUAL
temporal constraint gives

    W0=S0-sqrt(zeta)div pi_A/(a^3 sqrt(kappa)).
    ||Wbar||2 <=512delta U_psi+4e6delta integral_(t0)^t U_psi,
    ||S-Wbar||2 <=1024delta U_psi+4e6delta integral_(t0)^t U_psi.

Indeed the spatial and divergence energies together give
sqrt(2)*641*B/sqrt(zeta)<4e6 for the last time-integral coefficient.
This is why dropping J0 from the temporal solution is not allowed.

At a>=1 the physical inverse metric has Euclidean operator norm
at most1, while the spatial volume a^3 is at most B. Combine
this with ||D S[eta]||2<=2048delta U_eta. On the interval of
length1, Cauchy-Schwarz bounds the Volterra time integral as a
map L2 to L2 by1. Therefore

    |Force_psi[eta]|/kappa
      <= B*2048*(1024+4e6) delta^2
                       ||psi||H3(IxR3) ||eta||H3(IxR3)
      <4e10 delta^2 ||psi||H3 ||eta||H3.

For psi=epsilon f this is cubic in epsilon. The domain
delta<=1/100 alone does not make the normalized coefficient
small. A concrete smaller-ball corollary is

    delta<=1e-8  => |Force_psi[eta]|/kappa
                          <4e-6 ||psi||H3 ||eta||H3.

This is only the source/contact contribution in a stated dual
norm; it is not a bound on the full light inverse or a nonlinear
stability theorem. No relative error against the vanishing
classical bounce density is used.

Homogeneous sources give Wbar=S and cancel this force exactly.
The spatially nonclosed source control shows why such a
cancellation cannot be assumed for the whole perturbation class.
Neither case determines the metric/noise response of the
connected Gaussian determinant.
