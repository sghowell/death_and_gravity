# Physical readout matrices and local contact terms

Use the original physical canonical pair y=(w,p_w), not the
moving normalized oscillator pair. Write the mode Hamiltonian
as y^T K_H y/2, with K_H=diag(g² Omega_bare²,g^-2).

At fixed physical canonical data and comoving momentum,
the actual energy and pressure matrices are

    K_rho=(partial_N K_H)/a³,
    K_pressure=-(a partial_a-2q partial_q)K_H/(3N a³).

The N derivative includes both mass functions. The spatial
variation includes the physical momentum change q=k_com²/a².
These identities follow from the Legendre transformation
of the original metric variation.

On the clock, use the map to (v,p_w/g) to check them against

    rho_T: diag[q+m²(1+beta),1]/a³,
    rho_L: diag[(q+m²)(1+beta),1-alpha*z]/a³,
    p_T:   diag[(q-m²)/3,1/3]/a³,
    p_L:   diag[-(q+m²)/3,(1+2z)/3]/a³.

These are exactly the frozen physical readouts, not the
ordinary-Proca comparison energy in place of the actual source.

For a metric variation, the contact matrix is

    delta K=n partial_N K+zeta(a partial_a-2q partial_q)K.

The exact second lapse mass jets are derived from the full
frozen mass functions and p_affine(N), not an on-clock
constant-mass replacement:

    a_m,NN=-4/(3h)-8/(9h²),
    b_m,NN=-28/(27h)+152/(729h²).

The energy contact matrix depends on them. Suppressing that
matrix while varying only the mode covariance omits a real
part of the same action's second variation.

For symmetric equal-time canonical covariance Sigma,

    Sigma'=M Sigma+Sigma M^T,
    delta Sigma'=M delta Sigma+delta Sigma M^T
                    +delta M Sigma+Sigma delta M^T,
    delta O=Tr(delta K Sigma+K delta Sigma)/2.

The initial value delta Sigma(u0) is part of the chosen
state family. Both equations are checked as exact matrix
identities; neither is a numerical norm bound after UV
integration.

The retarded covariance and local contact contributions
must be combined with the corresponding second-variation
counterterms before assigning a finite response. S6.53
controlled first variations and specified their finite
prescription; it did not supply every required mass/curvature
second functional variation. The missing integrated
renormalization and function-space bounds remain work.
This follows the structural distinction discussed in
[S6.57's literature note](../../notes/literature.md); it does
not import a scalar flat-space kernel for this vector.
