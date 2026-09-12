# Literal inhomogeneous ADM scalar action

The regular physical/clock metric chart is S174's exact C=R^-1/2 map. Its full lapse-gradient and lapse-velocity cancellations include the boundary primitive I and its nonzero I_NN. S180's coefficient_clock_N_jets already contain that primitive. The present calculation uses those classical coefficients only; it does not use S180's old isolated quantum inverse.

At quadratic scalar order, rotational invariance permits a Fourier axis x. With v=v(t,x), contravariant shift beta=beta(t,x), define rate=H+v_t-beta v_x. The clock-chart extrinsic curvatures are
\[
 K^x_x=(rate-\beta_x)/N,\qquad K^y_y=K^z_z=rate/N,
\]
and
\[
 {}^{(3)}\hat R=a^{-2}e^{-2v}(-4v_{xx}-2v_x^2).
\]
The full transformed kinetic coefficient is R^(1/4)/2 multiplying KijKij-K^2. The linear trace coefficient is the full primitive-corrected b(N) from S180. The spatial curvature coefficient is R^(3/4)/2, whose first N jet is-3delta/2. The homogeneous potential is unchanged.

The original matter kinetic factor is U=R^-3/4. Its rate is ell+sigma_t-beta sigma_x; its quadratic spatial gradient is-sigma_x^2/(2a^2). The full canonical matter wave is retained at nonzero transfer.

Expand N=1+epsilon n, v=epsilon v, beta=epsilon beta and sigma=epsilon sigma through epsilon^2. The exact implementation differentiates the complete density, not a presumed reduced action. It then performs:

1. The weighted time integration by parts of the v v_t term.
2. The complete spatial divergence proportional to beta v_x+v beta_x.
3. The matter spatial boundary -ell beta sigma_x to+ell beta_x sigma.
4. The v v_xx spatial integration by parts, retaining n v_xx separately.
5. The weighted matter boundary3ell v sigma_t to-3ell v_t sigma, using ell'+3H ell=0.

For one real Fourier mode, v_x^2 maps to a^2q v^2, sigma_x^2 to a^2q sigma^2 and v_xx to-a^2qv. The resulting action per kappa a^3 is exactly
\[
 L_2=-3\dot v^2+(J+w^2/2-3\Theta^2)n^2+6\Theta n\dot v
 +\dot\sigma^2/2+wn\dot\sigma-3\ell\dot v\sigma
 +2b(\dot v-\Theta n)+\ell b\sigma
 +qv^2+2Eqnv-q\sigma^2/2.
\]
The full literal residual vanishes. The same J agrees with the independent S174 canonical lapse pivot. The shift equation is2(dot v-Theta n+ell sigma/2)=0 for each P nonzero.

An independent numeric fixture reconstructs the complete three-by-three extrinsic curvature, its trace and squared trace, full conformal curvature, matter transport and gradient. Its spatial mean agrees at both slab endpoints and the bounce, including momentum1e-9. The homogeneous part is subtracted only to isolate the spatial terms in that fixture; it is not removed from the actual action.

At P=0, the shift potential has no literal divergence variation. The nonzero-P equation cannot simply be imposed on the global homogeneous solution. No division by H, Theta or a vanishing background density is used.
