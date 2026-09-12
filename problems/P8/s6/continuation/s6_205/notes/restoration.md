# Restoring finite oversubtractions without changing matching

At the original finite regulator, partition the complete first-five-endpoint spatial Taylor sector exactly:

P_all,K = P_fin,K + P_UV,K.

Return the ten finite cells, with their actual signs, to the already controlled response:

J_actual,K = (known_S203,K + P_fin,K) + (C_unit,K + P_UV,K).

No finite term is discarded and no matching prescription is changed. The full actual-state/reference correction from S6.201 and the S6.203 all-momentum Taylor-subtracted remainder remain in known_S203.

Use

M[D]^2 = ||D||L2^2 + ||grad D||L2^2,
N61[Gamma]^2 = sum_(r=0)^6 ||partial_t^r Gamma||L2^2 + ||grad Gamma||L2^2,
X46[Gamma]^2 = sum_(r=0)^4 ||(1-Delta)^3 partial_t^r Gamma||L2^2,
Y[Gamma]^2 = N61[Gamma]^2 + X46[Gamma]^2.

The S6.203 known bound 4e54 M[D]Y[Gamma] and tail 6e60 M[D]Y[Gamma]/K, together with the finite-polynomial bounds, imply

|known_new| < 5e54 M[D]Y[Gamma],
|known_new - known_new,K| < 7e60 M[D]Y[Gamma]/K.

Both canonical metric factors supply 4/kappa, with kappa = 10^800. The finite polynomial alone has canonical bounds 4e-750 and 4e-745/K. The updated known actual piece has 2e-745 and 28e-740/K.

Every finite cell has positive spatial degree. Consequently this restoration is zero at exact P = 0 and preserves the earlier homogeneous result.

The fixed local spatial Hessian of S6.204 is a matching target, not an independently established quantum contribution. It is NOT added again here. Equality of the remaining quantum contact/Taylor sector to the original covariant prescription, including all finite coefficients and regulator artifacts, remains to be proved.
