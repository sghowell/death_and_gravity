# Two distinct momenta and the complete local contact

At gamma0, a detector Fourier component D_(k-q) gives

    M_D(k,q)=diag(
       -a m^2 D_(k-q)+C_k^t D_(k-q) C_q/a,
       D_(k-q)/a).

In particular the magnetic vertex contains BOTH C_k and C_q.
It is not generally symmetric as a single ordered pair; the
reversed pair has the proper transpose/adjoint. Real sources
satisfy the corresponding Fourier reality conditions.

For two metric directions D,Gamma, the exact mixed derivative
of exp(gamma) at zero is (D Gamma+Gamma D)/2. The same mixed
derivative of exp(-gamma) has the same positive sign.
Consequently the second contact has a positive mass term,
the same magnetic and electric form, and the full Fourier
convolution of this pointwise anticommutator. The constraint
term is independent of gamma and contributes no such vertex.

Define the ten-row energy feature map F_k on Z=(A,pi):

    F_k Z=(sqrt(a)m A, C_k A/sqrt(a),
           pi/sqrt(a), (k.pi)/(a sqrt(a)m)).

Then F_k^t F_k=M0(k) exactly. The first vertex is
F_k^t diag(-D,D,D,0) F_q and the contact is
F_k^t diag(H,H,H,0) F_q, H=(D Gamma+Gamma D)/2.
This proves the two-sided energy-relative bounds <=||D||op
and <=||D||op||Gamma||op for a fixed Fourier pair. For general
fields, the second bound is applied inside the full convolution.

For the actual symplectic normalizers
T_k^t M0(k) T_k=omega_k I, the matrices
F_k T_k/sqrt(omega_k) have orthonormal columns.
Therefore

    ||T_k^t M_D(k,q) T_q||op
        <=sqrt(omega_k omega_q)||D||op.

The proof is uniform at zero, noncollinear and arbitrarily large
momenta. It does not assume that normalizers at different
momenta commute. The extra constraint row is retained in the
energy identity even though its direct vertex is zero.
