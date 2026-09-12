# Complete weighted angular contraction

The removed leading shell has weight w(n) = (-n.phat)_+. Its required moments are

I0 = pi,
I2_ij = pi(delta_ij+phat_i phat_j)/4,

and

I4_ijkl = pi/24 [
  delta_ij delta_kl + delta_ik delta_jl + delta_il delta_jk
  + the six distinct delta_ab phat_c phat_d terms
  - phat_i phat_j phat_k phat_l ].

The code checks all 81 fourth-tensor components, not selected contractions. In an axis frame, each monomial moment follows from independent beta integrals in azimuth and u = -n3. Numerical quadrature checks every monomial of total degree at most four. Exact order validation occurs before caching.

Contracting the full leading Proca symbol gives

integral w(n) h(D,G,n) dOmega
 = pi [18tr(DG)-12(D phat).(G phat)
       -(phat.D.phat)(phat.G.phat)]/(64a).

The Fourier measure (2pi)^-3 changes the denominator to 512 pi^2 a. No extra Wick, retarded or metric factor is inserted: the current normalization is the one already fixed by the complete endpoint formula.

In the Frobenius-orthonormal basis of two tensor, two vector and one scalar tracefree spatial matrices, the bracket is diagonal with eigenvalues 18, 18, 12, 12 and 28/3. The code verifies the full 5-by-5 matrix and orthonormality, and independent tests rotate this entire basis. Hence the coefficient is positive on every nonzero real tracefree D = G, not only on a TT slice.

The LL contribution changes this coefficient. Omitting it is a different two-polarization calculation, explicitly rejected as a control.
