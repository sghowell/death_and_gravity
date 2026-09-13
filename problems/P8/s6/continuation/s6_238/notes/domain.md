# Full real domain, regular affine maps and the upper edge

The original R=1+B(X)(X-1)/(1+u²)^3 obeys L<R<6/5, L=1/2+1/(8N), N=1024, for every real u on -1/4096<X<6/5. Here

T=X^N/[X^N+(1-X)^N],
bump=N X² exp(-N X²),
B=T+(1-T)bump.

Because N is even, 0<=T<=1 and0<=bump<1. The new addition is deltaR=N X²(1-X)^8 exp(-A X²)>=0. On the whole original strip, |1-X|<2 and X²exp(-A X²)<=1/A, so deltaR<=2^18/A<10^-414. The lower bound is unchanged.

A tiny additive bound alone does NOT preserve a strict upper bound without a uniform margin. For X<=1, R_old<=1; for1<=X<=9/8, R_old<=X<=9/8. These ranges plus2^18/A stay strictly below6/5.

For9/8<=X<=6/5, one has1-bump>1/2 (the maximum of v exp(-v) is1/e<1/2), and

1-T>=(1/2)(5/48)^N,
(X-1)(1-B)>=(5/48)^N/32.

Since (1+u²)^-3<=1, R_old<=X-(X-1)(1-B)<=6/5-edge_gap. Meanwhile deltaR<=2^18 exp(-A)<=2^18*4!/A4<edge_gap/2 by the positive exponential series. Thus the SAME upper6/5 persists uniformly, even at the endpoint limit. No old-domain theorem is applied outside its range.

The full R_new-1 has an X² factor with analytic coefficient. Therefore the generic S174 covariant chart

g_phys=R^-1/2 g_hat+(1-R^-1/2)du du/X

has its regular divided coefficient, unchanged scalar norm X_phys=X_hat, determinant ratio R^-3/2 and nonzero field-space Jacobian R^-9/2. The inverse divided coefficient is regular as well, including nonzero null gradients. The physical matter metric and cone are not changed by fiat.

Use the NEW unique vacuum-regular lower solution

q=(3X/4)R(u,X)^(3/4) integral_0^1 sqrt(t) R_X(u,tX)R_u(u,tX)R(u,tX)^(-7/4)dt.

The segment tX stays in the connected strip and R remains positive. This parameter integral is analytic on compact subdomains, from either sign of X. If R=1+a(u)X²+..., then q=a a_u X4/3+..., solving
q_X+[1/(2X)-3R_X/(4R)]q=3R_XR_u/(4R).
Its clock value is generally not zero and must not be reset.

The already proved GENERIC full rank60 connection identities in S174 apply to this new R range. With p=sqrt(R)/2, 1/3<p<3/5, tau=3(2p³-1)/p<0 and sigma=(8p+5)/(8p²)>0. Recompute Xi=eta-D_old^-1 and retain the entire source-centering terms. The exact full quotient determinant ratio is -tau sigma³>0; the56 complement directions remain algebraic and the four projective directions gauge. The original generic all60 Euler-lift and all64 source identities are replayed as frozen inputs, not replaced by a vector-only ansatz.

Recompute B_shift=3Href(R-1)+3(q+R_u/2)/2 and W=T_trace-B_shift du. The q dependence cancels identically from the retained source

S_mu=u_mu(R-1)[-3Href+3R_u/(4R)+Box u/X+(-1/X²+3R_X/(2RX))uHu].

The new q, B_shift and affine background trace need not equal the old ones. Add the heavy scalar and new f independently of the connection. The full generic construction then yields the desired NEW reduced covariant action. It does not prove a healthy physical constraint rank or characteristic cone on every field background; those are separate dynamical questions.
