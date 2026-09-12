# Formulation and exact boundary

Keep the same S6.176 actual CD canonical Proca sector, m=1000, kappa=10^800, fixed prepared state, unit time slab and original finite renormalization prescription. The complete actual-to-unit-W8 comparison correction is retained from S6.201. No state or prescription is changed.

Let Amax=25/16, nu(k)=sqrt(m^2+|k|^2/Amax^2), delta=10^-6 and E=2e27. Write E5(P,k) for the complete PRE-CURRENT row of the first five unit-W8 upper endpoints, from source time jets0,...,4 to the full detector tensor. The final imaginary part is applied only after the row is assembled. All nine physical pairs and the constrained temporal readout remain.

At the original two-leg regulator chi_K=1_(|k|<=K,|-k+P|<=K), define

Q_K = integral chi_K [1_(|k|<m) E5(P,k)
                      +1_(|k|>=m)(E5(P,k)-T4 E5(P,k))] d^3k/(2pi)^3,

with test contractions, Fourier integration and time integration understood. T4 is the degree-four external-spatial Taylor polynomial of the coefficient at fixed k, not of the tests or regulator. This fixed internal split is a mathematical comparison partition, not a physical cutoff.

The high band is partitioned for estimation into FAR: |P|<=delta nu/2 and NEAR: |P|>delta nu/2. S6.202 supplies the far bound. Here the full real endpoint row, unexpanded low band and complementary near region give

|Q_K|,|Q| <3e54||D||L2 X46[Gamma],
|Q-Q_K| <5e60||D||L2 X46[Gamma]/K,

for K>=m, where X46^2=sum_(j=0)^4||(1-Delta)^3 partial_t^j Gamma||L2^2. These are absolute weak bounds and an explicit convergence rate with the original removed union.

Let N61^2=sum_(j=0)^6||partial_t^j Gamma||L2^2+||grad_x Gamma||L2^2, Y^2=N61^2+X46^2, and M[D]^2=||D||L2^2+||grad_x D||L2^2. Combining S6.201 gives the known actual-current finite piece R_actual-unit+F_unit+Q below4e54 M[D]Y[Gamma], with error6e60 M[D]Y[Gamma]/K. Its canonical bounds are16e-746 and24e-740/K.

The exact remaining sector is C_unit,K+P_K, where P_K integrates the high-band Taylor INTEGRAND with the original two-leg chi_K and C_unit,K keeps its distinct one-mode band. P_K need not be a polynomial of P after the moving band is integrated. No finite/divergent coefficient, regulator artifact or fixed covariant matching of this sector is assumed away.

No full matched response, derivative-compatible mixed inverse, interacting background, stability, physical cutoff, remaining parent loops or finite-gravity IR/Regge result follows. Original V/G/B and P8 remain OPEN.
