# Full complex-disk coefficient tails and the light logarithm

The estimates apply for any fixed real |L|,|b|<=K. On the complex circle |w|=1/(16K),

abs(-4Lw+4Lb w²)<=17/64.

Choose the analytic square root d(0)=1. Its binomial series, or the exact square-root difference with its analytic positive-real branch, gives |d-1|<1/4. Indeed 17/64<7/16=1-(3/4)². Thus |d|>3/4. For a=(1-2Lw+d)/2,

|a-1|<1/16+1/8=3/16, |a|>13/16,
|Lw/a|<1/13<1/8,
|beta|<2/13<1/4.

The argument (1+d)/2 lies within1/8 of1. Its analytic logarithm is bounded by1/7<1/4. The entire H series gives

|H(beta)|<=1+|beta|/[2(1-|beta|)]<2.

Therefore the complete analytic coefficient functions in the exact primitive obey |P|<2 and |Q|<4 on the full circle and inside it. In particular |P|<3/2 and the conservative Q estimate is73/24<4. These are continuum bounds; numerical circle samples are only independent diagnostics.

Let R=1/(16K), u=w/R<=1/2 and keep the explicit H0=-ln(w)-Log(L-i0). Cauchy's coefficient bounds give, for the complete tails after degree1,

|P_tail|<=4u², |Q_tail|<=8u².

The exact weighted geometric sum is

sum_{k>=2}(k+1)u^k=u²(3-2u)/(1-u)²<=8u².

It controls derivatives of ALL the tail coefficients. Hence

|P_tail+w P_tail'|<=16u²,
|Q_tail+w Q_tail'|<=32u².

Using -partial_n=w² partial_w and partial_w H0=-1/w, the full remainders satisfy

|R_T|<=4w u²(|H0|+2),
|R_D|<=w²u²(16|H0|+36)<=16w²u²(|H0|+3).

This differentiation includes the explicit logarithm; Cauchy is never applied to -ln w at its singular origin.

## Uniform physical and symmetric parameter domains

For a real physical configuration with S>=4, the channels are S and two numbers in[-(S-4),0]. For every xi,eta in[0,1],

|1-channel*xi(1-xi)|<=S/4,
|other_channel*eta(1-eta)|<=S/4.

Thus take K=S/4. If S/n<=1/8, then u=4S/n<=1/2 and

|R_T|<=64S²n^-3(|H0|+2),
|R_D|<=256S²n^-4(|H0|+3).

At the symmetric subthreshold point take S=4,K=1; its actual channels4/3 obey the same bounds.

## Integrable threshold and full physical boundary

For the physical s=S>=4, factor L=S(xi-r_plus)(xi-r_minus), where both light roots lie in[0,1]. For any r in[0,1],

integral_0^1 -ln|xi-r| dxi
 =1-r ln r-(1-r)ln(1-r)<=1+ln2<2,

with endpoint x ln x limits included. Therefore

integral |Log(L-i0)|<=ln S+4+pi.

This includes the coincident roots at S=4 and every intermediate physical energy. The imaginary contribution is bounded by pi; its sign is retained rather than discarded.

For either negative channel, L is between1 and S/4, so the logarithm is simpler and satisfies the same majorant. At the symmetric point,2/3<=L<=1 and |ln L|<1. On the actual model n<10^198 and ln10<7/3, hence ln n<462. The latter elementary bound follows already from the positive exponential partial sum through degree8 at7/3 being larger than10.

Since S<n, both integral |H0|+3 and |B(s)| are below1000, uniformly over all physical angles and the separate symmetric point. The resulting integrable logarithmic majorant justifies the parameter boundary limits, all integrations of the coefficient tails, and the displayed complete box derivative. No ordinary absolute-value bound is applied to the original Feynman double-pole integrand.
