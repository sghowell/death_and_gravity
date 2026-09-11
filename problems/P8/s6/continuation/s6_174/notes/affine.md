# Exact full-connection mass change and retained vector

Use the S6.109 analytic affine inverse and its unique regular q.
All sources refer to the unrestricted connection, with four projective
gauge directions and a rank60 quotient. Let y be the displacement
from its actual stationary connection, M the quotient Hessian, and
N the four-trace map for T=V_trace+U_trace. Its full inverse response is

    D_old=N M^(-1)N^T=diag(tau_p,sigma_p,sigma_p,sigma_p),
    tau_p=3(2p^3-1)/p, sigma_p=(8p+5)/(8p^2), p=sqrt(R)/2.

The displayed frame is source-signature-+++ orthonormal; use
eta=diag(1,-1,-1,-1), the P8 physical Lorentz matrix.
On1/2<R<6/5, p is between1/3 and3/5, tau_p<0 and sigma_p>0.
Neither D_old nor the already certified full quotient is singular.

Choose the NEW tensor Xi=eta-D_old^(-1), not constant(11/9)eta.
Add one half(T-Tstar)^T Xi(T-Tstar), retaining its linear and
constant source-centering terms. It equals the old mass change at
R=1 but differs away from that surface.

Covariantly in P8 signature Xi=a g^(-1)+b du^sharp du^sharp,
where a=1+1/sigma_p and b=(-1/tau_p-1/sigma_p)/X.
The numerator of b has the factor2p-1=(R-1)/(2p+1).
Thus b=O(X) analytically at X=0, including nonzero null gradients.
There is no preferred-vector direction singularity hidden there.
The source-pinned affine coefficients and q are analytic on the
same connected domain. The timelike tensor identities extend
analytically to the full domain; they are not a rest-frame truncation.

Put L=M^(-1)N^T D_old^(-1). Then N L=1, and the full updated
matrix obeys

    M_new=M+N^T Xi N,
    M_new L eta=N^T, N L eta=eta,
    L^T M_new L=eta.

The projector1-LN annihilates N and has zero mixed quadratic
term with L. Its56-dimensional kernel sector is unchanged and
nondegenerate. Therefore its elimination is exact, leaving
one half(T-Tstar)^T eta(T-Tstar). The determinant ratio is
det(1+Xi D_old)=det(eta D_old)=-tau_p sigma_p^3>0; no new quotient
zero is introduced. The four projective directions remain gauge.
The code checks all60 Euler-lift rows, not only a four-vector ansatz.

With the specified curl of W=T-B du, the exact retained action is

    S_target+S_M1+kappa int sqrt(-g)[
        -zeta F(W)^2/4+(W-S)_mu(W-S)^mu/2].

All56 complementary connection directions have been eliminated
without a derivative expansion. W remains dynamical; it has NOT
been integrated out. Its nonlinear source S is specified separately.
The algebraic identity does not establish the quantum measure,
loop stability of Xi or a UV completion.
