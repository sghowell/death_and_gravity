# Full shift transport and its sharp complex Fourier norm

The reduced shift term is

    integral beta^i[Pi^j F_ij-A_i divPi].

For a covector A,

    (Lie_beta A)_j=beta^i partial_i A_j+A_i partial_j beta^i.

Their difference is the explicit divergence partial_j(beta^i A_i Pi^j). Thus the complete shift Hamiltonian generates covector Lie transport modulo the compact spatial boundary. The longitudinal term-beta.A divPi is essential; retaining only the familiar magnetic-electric flux would not generate this transport.

At P=k-q, the complete Fourier blocks in Z=(A,Pi) order are

    M_PiA=i[(beta.q)I+P beta^t],
    M_APi=-i[(beta.k)I-beta P^t].

At P0 both canonical field components translate with i beta.k. The formula retains distinct k,q and the reverse-pair adjoint condition.

In the ten-field ordering of vertices.md, the only nonzero shift blocks are

    S_BE=S_EB=-i a C_beta,
    S_AL=-i a beta, S_LA=+i a beta^t.

The magnetic/electric and mass/divergence blocks act on disjoint groups. Multiplying F_k^t S_beta F_q gives precisely the full Lie blocks above; it avoids a spurious separate |P|/m bound.

For real beta, S_beta is Hermitian, S_beta^3=a^2|beta|^2 S_beta and tr(S_beta^2)=6a^2|beta|^2. But Fourier beta may be complex, so that real spectral identity alone would be insufficient. For arbitrary complex beta set T=S_beta^dagger S_beta. Direct complex cross-product algebra gives

    T^2=a^2||beta||_2^2 T,
    trT=6a^2||beta||_2^2.

T is positive Hermitian, hence its nonzero eigenvalues all equal a^2||beta||_2^2. For beta!=0 its trace is positive, proving the sharp operator norm

    ||S_beta||op=a||beta||_2,

with zero at beta0. There is no additional sqrt2 or external-momentum factor. Exact six-real-component identities and independent complex-direction singular values check the statement.

In coordinates(N,beta^i,Q) the reduced Hamiltonian is linear in beta and independent of N,Q in its shift part. Therefore its shift-shift and shift-other second Hamiltonian vertices are zero. This is not deletion of covariant metric contacts: the map to g_mu_nu is nonlinear and has nonzero second derivatives.
