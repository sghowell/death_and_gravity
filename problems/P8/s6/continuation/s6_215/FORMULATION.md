# P8 S6.215: prepared Ward reconstruction and three ordered scalar kernels

This checkpoint reduces the remaining full conditional Gaussian metric response problem. It gives the exact ordered Ward reconstruction, correct source/detector endpoint conditions, full one-current chart correction, and explicit derivative-losing bounds on the known pieces. It does not construct the remaining scalar kernels.

The parent, physical signature, mass1000, kappa1e800, all-order initial state and S176 covariant finite prescription are unchanged. Sources vanish on the common initial neighborhood; detectors vanish on a final neighborhood. Begin with smooth compact variations and spatial Schwartz data. The initial-only source extension needed by the synchronous map is proved in notes/boundary.md.

For g=diag(1,-a^2 I), a=(1+t^2)^2, split X=Xsyn+Xxi with n_xi=eta', beta_xi=chi'-a^-2 grad eta, Q_xi=2H eta I+symgrad chi. Use initial retarded primitives for G and final advanced primitives for D. The full response is its synchronous spatial part plus two distinct ordered Ward terms involving the actual one-point density and complete S214 nonlinear metric chart.

Rotation covariance leaves exactly three additional ordered scalar kernels beyond S213. Their existence, original-prescription matching and quantitative bounds are not supplied here. The conditional reconstruction estimate is not a full response theorem, reduced inverse or quantum background result. Original V/G/B and P8 remain OPEN.
