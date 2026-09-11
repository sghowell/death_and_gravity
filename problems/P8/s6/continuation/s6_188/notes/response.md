# Covariance tangent, Kubo relation and metric contact

Write Z=(A,pi), J=[[0,I],[-I,0]] and A_gamma=J M_gamma.
For the symmetrically ordered real covariance
C=Re<ZZ^T>, the exact finite-mode equation is

    C'=A_gamma C+C A_gamma^T.

For a prescribed perturbation gamma=epsilon eta(t)D
that is zero near the common Cauchy surface, its
derivative at epsilon0 obeys

    delta C'=A0 delta C+delta C A0^T
       +eta[J M_D Cbar-Cbar M_D J],
    delta C(t0)=0,

where Cbar(t) is the actual unperturbed state covariance.
No quantum metric initial-state factorization is assumed.

Varying the coordinate current J_G=-<H_G> gives

    delta J_G=-1/2 tr(M_G delta C)
               -eta/2 tr(M_GD Cbar).

The second term is the metric-variation contact. It is
not removed by centering a fixed-background stress.

In the same finite-mode algebra, [Z_a,Z_b]=iJ_ab and
H_M=Z^T M Z/2 imply

    [H_M,H_N]=(i/2)Z^T(M J N-N J M)Z.

Transporting the readout to the source time and using
trace cyclicity makes the covariance tangent identical to

    delta J_G(t)=-<H_GD(t)>eta(t)
       +i integral_(t0)^t <[H_G(t),H_D(s)]>eta(s) ds.

The sign follows from perturbing the Hamiltonian by
+epsilon eta H_D and reading out -H_G. Independently,
with W=C+iJ/2, Wick's theorem gives
<H_M H_N>connected=tr(M W N W^T)/2; subtracting the
two orderings gives the same commutator and factor.

Primary-source context: the distinct commutator and
second-metric-variation kernels in
[Hu and Verdaguer,0802.0658v1,eq.(4.6)](https://arxiv.org/pdf/0802.0658)
support keeping both contributions. Their scalar kernel
is not substituted for the Proca calculation above.

These statements hold modewise before the UV limit.
Integrating them into the fixed covariantly renormalized
stress response still requires its subtractions and
local finite metric terms. The present contact is not
a claim to have computed all renormalization contacts.
