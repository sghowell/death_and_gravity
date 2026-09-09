# Exact finite-frequency canonical system

Let x=(b,chi), with old momentum P in the S6.89 high-momentum
coordinate chart. Its Hamiltonian Hessian blocks are A,E,C, with
x'=AP+Ex. Put B=alpha/r_N. The canonical shift

    P=Pi+S_q x, S_q=q diag(B,0)

retains the weighted contact q diag(B'+Hhat B,0).
Directly expanding the literal full Hamiltonian, finite.py verifies

    A=A0+A1/q+A2/q^2,
    Enew=E+A S_q=E0+E1/q,
    Cnew=C+E^T S_q+S_q E+S_q A S_q
         +q diag(B'+Hhat B,0)=qG+C0.

A0 is exactly the inverse of S6.89's principal K and G is its
Euler-first principal gradient. Every remainder matrix is retained.
This construction uses no finite-q Legendre inverse: the old chart's
possible Legendre-Hessian zeros do not obstruct this canonical
first-order evolution. The only chart restrictions here are
q>0 and r_N!=0, already satisfied in the packet band.

Normalize Q=R^(3/2)x and P_c=R^(3/2)Pi to remove the weighted-volume
canonical damping. For k=|k_vector| constant in comoving coordinates,
define Y=(kQ,P_c). The exact Fourier ODE is

    Y'=[kJ+L0+L1/k+L2/k^2+L3/k^3]Y,
    J=[[0,A0],[-G/R^2,0]],
    L0=diag(E0+3Hhat/2 I,-E0^T-3Hhat/2 I),
    L1=[[0,R^2 A1],[-C0,0]],
    L2=diag(R^2 E1,-R^2 E1^T),
    L3=[[0,R^4 A2],[0,0]].

All blocks are exact and have the constant symplectic structure.
The positive symmetrizer diag(G/R^2,A0) cancels the leading kJ
in its energy derivative. There is no exponential-of-k estimate.

## Action-normalized modes

Use v_c=(1,-ell), v_m=(0,1), with kinetic diagonals
kappa_c=-h/(4m_N^2), kappa_m=e^3/N. Their positive coordinate
frequencies are omega_c=N c_clock/(eR), omega_m=N/(eR).
Columns of the mode matrix are

    S_j^sign=(v_j, sign*i*omega_j*K*v_j)
             /sqrt(2omega_j*kappa_j).

The order is (+c,+m,-c,-m). The eigenvalue matrix is i Lambda,
Lambda=diag(omega_c,omega_m,-omega_c,-omega_m).
For the real-time symplectic form Omega,

    i S^dagger Omega S=diag(-1,-1,1,1).

modes.py verifies this and the complete eigenbasis directly for
independent kinetic/frequency parameters. With Y=SU, the leading
transport is

    B0=S^-1 L0 S-S^-1 S'.

The calculation retains derivatives of ell, both kinetic diagonals
and both frequencies. All four diagonal entries of B0 vanish
identically, and its Krein symmetry is also checked. This relies
on the full real block-diagonal L0 and the action normalization;
it is not imposed by deleting amplitude terms. Consequently the
large crude norm of B0 will not become a spurious diagonal
exponential in the wavepacket error proof.
