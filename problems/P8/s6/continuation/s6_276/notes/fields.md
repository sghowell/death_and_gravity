# Entire canonical field realization and connected input path

Letx,y,z be original torus coordinates, k=P=10^64, and0<=s<=s*. Set

    v=sqrt(6s)/k cos(kx), Q=I, tau=0, Pi_v=0,
    Pi_tau=kappa sqrt(s/2) cos(ky)diag(1,0,-1),
    M1 spatially constant,
    Pi_M/kappa=1/10+d cos(kz).

Here d is the small fixed number derived in notes/fold.md. All heavy and
vector coordinate/momentum inputs vanish. These use only the original
axis sine/cosine modes and retain the original homogeneous means.
The scalar and tensor input amplitudes vanish continuously at s0.

The S267 base cotangent lift gives the EXACT reconstructed fields

    gamma=e^(2v)I,
    pi=-kappa e^(-2v) sqrt(s/2) cos(ky)diag(1,0,-1).

Indeed DQ* applied to Pi_tau gives precisely this pi, and2pi:gamma=Pi_v=0.
The exact density-weight-one momentum generator is

    D_i=-2partial_j(gamma_ik pi^jk)+pi^jk partial_i gamma_jk
         +Pi_M partial_i M1 + full vector/heavy terms.

The tensor matrix is tracefree and transverse to its y-frequency.
The first metric term is its ordinary divergence and the second is
proportional to its trace. Both vanish; all matter terms vanish because
their coordinates are constant or zero. Thus D0=0 identically. The full
mean-zero adjoint correction is exactly zero, not neglected. SinceQ=I,
the spatial gauge operator is the flat mean-zero invertible operator,
independent of the nonconstant conformal factor. No global quotient
theorem is inferred. Every residual translation charge
integral(Pi_v gradv+Pi_tau:gradtau+Pi_M gradM1+...) is zero as well.

The full Christoffel/Ricci construction, checked independently of the
conformal formula, gives

    c=e^(-2v)[4k^2 V cos(kx)-2k^2 V^2 sin(kx)^2],
    sh=2A^2 e^(-6v)cos(ky)^2,

whereV=sqrt(6s)/k andA=sqrt(s/2). At x=pi/(2k),y=z=0 one hasv0,
gammaI, c=-12s, sh=s and pm=sqrt(m*^2). The other original invariants
are zero. All three cyclic rotations of this construction are checked.

At s0 the metric isI and the only nonzero deviation is the M1 momentum
d cos(kz). Every invariant at EVERY point lies in the old S266 box,
so the original local positive-lapse branch exists on the initial whole
spatial datum. At the designated point it is exactlyN1.

For larger s no claim of a full lapse field or time solution is made.
Curvature need not be small even though v is small at this high frequency.
The reconstructed exponentials contain every generated harmonic. The
proof does not project these fields back to the input Fourier band.

These are canonical test inputs, not simultaneous quantum eigenvalues,
an interacting state or a path traced by the physical time evolution.
The fixed-reference boundary shear of S275 is invertible and independent
of N, so it merely relabels these finite phase inputs; it cannot remove
a zero in C_N or annihilate all transverse derivatives of C.
