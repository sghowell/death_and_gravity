# Complete current and local Euler variation

At any evaluation time a constant canonical GL transformation can set
K0=I. Its matrix is fixed while taking time derivatives, so K1..K4 remain
arbitrary symmetric, generally noncommuting jets. The scalar frequency
jets omega0..omega4 are independent, with omega0>0. Tensor covariance
under this constant transformation transfers an identity at K0=I to
every positive K0.

At K0=I the balanced Hamiltonian variation is

    deltaG=diag(2deltaomega/omega I-deltaK,deltaK).
    J=-deltaomega tr(Sigma_QQ)
        +(omega/2)tr[deltaK(Sigma_QQ-Sigma_PP)].

Construct B time jets directly from B^2=K, with B0=I, retaining every
ordered product. Construct the inverse, L,R,S and finite Riccati
r1..r4 similarly. Finally form the complete covariance coefficients
from C=(I-rsharp r)^-1 and the QQ/PP blocks. Coefficient adjoints do
not conjugate the formal marker. This is precisely the same
fourth-order physical current comparison as in S6.192.

Now vary the root-free action independently. Set
A=K'K^-1, p=omega'/omega, s=pI/2-A/2 and t=s'-p s.
In particular deltaA=deltaK' K^-1-A deltaK K^-1.
Compact integration by parts gives

    E_K2=K^-1 t/(4omega),
    E_omega2=-tr(s^2+t)/(4omega^2).

For fourth order define U=-t'+2p t+2s^3,
v=trU-2tr(ts) and Z=tr(t^2+s^4). Full variation gives

    E_K4=Sym{K^-1[U'+[U,A]-3pU]}/(16omega^3),
    E_omega4=-[3Z+v'-3p v]/(16omega^4).

Symmetrization is required because deltaK is symmetric. The commutator
must not be discarded. These are raw local Euler derivatives, not
derivatives taken after holding the time jets fixed incorrectly.

The code equates EVERY word coefficient of E_K with
omega(Sigma_QQ-Sigma_PP)/2 at orders0,2,4, and every independent trace
coefficient of E_omega with -trSigma_QQ. There are4 matrix words at
order2 and16 at order4. It separately checks all reference symmetry
and real-phase word coefficients. Root/inverse checks require their
entire residual word maps to be empty, not just a sum of coefficients.

This proves the compact local variational identity for arbitrary
noncommuting positive K and independent frequency, hence for every
admitted physical metric history and detector by the chain rule.
Independent literal local-action Euler derivatives with both unit
and nonunit K0 supplement the exact word calculation.
