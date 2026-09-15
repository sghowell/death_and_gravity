# Actual reference symmetry and all homogeneous equations

This argument concerns the original prepared state, not an arbitrary
translation-invariant state. Translation invariance by itself does not justify
deleting directional, internal or imaginary cross covariance.

## The full finite representation

Let G be the24 proper signed-permutation matrices in three dimensions. The
original nonzero shell consists of the cosine and sine of P x_i, i=1,2,3.
For each axis and real wave the eight configuration channels are
(v,M1,tau_plus,tau_cross,H,W_1,W_2,W_3). For g e_i=epsilon e_j, cosine has
sign1 and sine has sign epsilon. Scalars transform by1, the vector by g,
and tensor coefficients by
T_ab=Tr[E_{j,a} g E_{i,b} g^T], with both orthonormal transverse traceless
matrices E retained. This gives U_g on48 configurations and diag(U_g,U_g)
on96 canonical phase coordinates. symmetry.py evaluates every signed action,
all576 compositions, orthogonality and the full symplectic identity.

The covariance used in that identity retains the arbitrary symmetric scalar
4x4 block on(v,M1,p_v,p_M1), a complete2x2 block for each identical tensor
polarization, a complete heavy2x2 block, and different complete transverse
and longitudinal Proca2x2 blocks. Both transverse vector blocks agree; the
longitudinal block need not equal either. Internal and q-p entries are not
set to zero. The full Proca matrices
K=a^-1 I+k k^T/(a^3 m^2) and
V=a m^2 I+(|k|^2 I-k k^T)/a
are explicitly covariant under all24 rotations.

## Binding to the unchanged seed

The actual S251 phase.py canonical generator, positive selection form and
both time-dependent canonical boundary maps depend on momentum only through
q=|k|^2/a^2. They are identical for all three axes and cosine/sine members.
The S251 gaussian.py minimizer of the full positive preparation Gramian is
unique, including noncommuting internal and q-p blocks. Thus the actual
selected scalar covariance has exactly the radial block structure above.
The two tensor polarizations have the same full generator and preparation.

The S240 state.py scalar preparation is an exact basis-independent SLE
depending on the radial frequency and the unchanged time bump. Its full
heavy covariance is the same for all directions and both real waves.
The S190 initial.py data retain the original all-order T/T/L Cauchy
frequencies and slopes. On the actual zero-shear reference these are radial,
with equal transverse blocks. The full frame equations and uniqueness of
their exact mode/Riccati evolution preserve that structure. No finite WKB
replacement, state average or isotropic reset is made here.

Consequently the actual physical seed covariance V0 obeys
S_g V0 S_g^T=V0. For its existing whitening V0=S0 S0^T/2, the action
S0^-1 S_g S0 is both symplectic and orthogonal. Hence the same radial
phase cutoffs, the full Laplacian D=Delta/4 and the same-state coherent
quantization are equivariant. Weyl covariance gives the other ordering.
This does not require changing S0 or discarding its cross covariance.

## Equivariance of reconstruction and evolution

Each full original scalar contraction, spatial derivative, Fourier matrix
projector, determinant/shape fixed point, mean-zero formal-transpose inverse,
temporal root and lapse root is equivariant. Their domains are invariant;
uniqueness of their local solutions preserves equivariance. Spatial
integration preserves it as well. Generated harmonics and reconstructed
means remain in the complete nonlinear symbols. The fixed free reference
flow is radial and commutes with these actions. Translations are also
preserved, including the original zero-total-translation-charge sector.

The centered pure Gaussian wavefunction is invariant under the ordinary
physical configuration pullback: covariance invariance determines it up to
phase, and its nonzero value at the origin fixes that phase to1. The fixed
whitening merely conjugates this representation. Both bounded full
Hamiltonians therefore preserve the invariant Hilbert subspace.

To discuss omitted homogeneous directions, first extend the same local
absolute canonical chart in the three spatial-vector and five traceless
symmetric shape coordinates and their conjugate momenta. At zero the full
auxiliary pivots and mean-zero inverse are invertible; their strict margins
give a smooth equivariant extension in these finite additional directions.
No quantitative excursion in those directions is asserted or required.
At invariant scalar homogeneous data and invariant quantum state, every
derivative of the energy in a vector direction is an invariant vector;
every shape derivative is an invariant traceless symmetric matrix. The
exact group averages satisfy
sum_g g/24=0 and sum_g g A g^T/24=Tr(A) I/3.
Thus all three vector and all FIVE shape coordinate AND momentum equations
vanish at zero. This is not a two-TT-zero-mode truncation. The heavy
homogeneous scalar is not removed by this symmetry. The temporal normal
vector remains the complete algebraically reconstructed source.
