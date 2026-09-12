# Exact analytic normalization at inverse radius zero

Write x = 1/r and use scaled momenta v_k = n and v_l = -n+xP. For either leg,

Omega = sqrt(v.v/a^2+m^2 x^2),
z = 1-m^2 x^2/Omega^2,
What = Omega + sum_(q=1)^4 P_q(t,z) x^(2q) Omega^(1-2q).

The normalized amplitudes are

F = 1/sqrt(2What),
Pi = [-i What - x(What_time/(2What)+rate)] F,

with the opposite i sign for the analytic Schwarz detector. The total time derivative retains the derivatives of Omega, z and all four actual W8 coefficients.

At x = 0, the two scaled momenta are n and -n, both projector denominators are one, and Omega = What = 1/a. The actual positive scale factor on the compact slab ensures that all roots and inverses have their chosen analytic branches in a sufficiently small complex disc. Compactness in time and the unit angular sphere makes the disc uniform for every fixed bounded external-momentum set. This does not claim a single momentum-independent radius or an explicit all-P norm.

Substituting f = sqrt(x)F, p = Pi/sqrt(x), omega = Omega/x and the scaled momenta into the exact sector expressions gives

x^2 F_pair,original
 = F_pair,normalized(n,-n+xP; algebraic mass m*x).

The code checks all four complete two-time identities. The effective m*x is only a normalization device: the original physical mass remains 1000.

All dependence on external momentum is through xP, and all inverse-radius denominators are regular at zero. No global smooth individual polarization frame is used. No analyticity of the prepared Borel state is assumed; the construction concerns only the unit-W8 comparison and retains the complete actual-state correction.
