# Canonical phase columns, not frozen independent oscillators

Let kappa=sqrt(q_initial) and a_i be the original global
scale factor at the left endpoint. Write the actual
positive kinetic matrix as alpha_i=B^T B with real
invertible upper-triangular B. With S_i the symmetric
part of beta_i, take the two initial column matrices

Q_i=B^-1/sqrt(2*kappa*a_i^3),
P_i=-i*sqrt(kappa*a_i^3/2)*B^T+a_i^3*S_i*Q_i.

They obey Q_i^dagger P_i-P_i^dagger Q_i=-i I and
Q_i^T P_i-P_i^T Q_i=0. The equal-time field and momentum
commutators also have their canonical values. The exact
matrix checks keep a general real B and symmetric S_i.
Real Hamiltonian evolution preserves these identities.

Crucially the initial velocity is

Q'_i=-i*sqrt(kappa/(2*a_i^3))*B^-1
     -alpha_i^-1*A_i*Q_i.

The antisymmetric second term is retained. Shifting P_i
by the full beta_i instead of S_i would remove this term
but spoil canonical isotropy. With alpha=I, kappa=a_i=1,
S=0 and A=[[0,1],[-1,0]], the isotropy defect is exactly A,
not zero. This is a tested negative control.

These finite-band data specify canonical mode columns;
they do not assert an all-momentum scalar Hadamard state
or an infinite-volume implementation theorem.

## Concrete energy-to-column inequalities

The original global scale satisfies 1<=a<=25/16 and
a^3<8 on I. From alpha>=I/1000, ||alpha||<10^4,
V<=2*10^4*q and kappa>=1000*||A_i||,

||Q_i||<32/sqrt(kappa),
||Q'_i||<64*sqrt(kappa),
E_i<10^9*kappa.

The preceding exact energy bound gives E<2*10^9*kappa
on the declared window. Scale drift makes q>=kappa^2/2,
so E>=||Q'||^2/2000+q*||Q||^2/4000 yields

||Q||<4*10^6/sqrt(kappa),
||Q'||<2*10^6*sqrt(kappa).

In constant spatial coordinates normalized at the
window center, take external magnitudes between
L0=10^12 and U0=10^13 and hard proper subset transfers
at least L0. Internal tree momenta are at most 2U0.
For all these modes, both initially and throughout the
window, L0/2<=sqrt(q)<=4U0. This overestimates the actual
scale drift. It also exceeds both chart q_min thresholds.
Consequently

||Q||<8, ||Q'||<10^7*U0,
||P_chart/a^3||=||alpha Q'+beta Q||<10^12*U0^2.

Here ||beta||<10^5*(16*U0^2+1); both its leading q term
and its remainder are included. The exact rational
margins for these displayed relaxations are recorded.

## Return to the full reduced Hamiltonian variables

For the unitary chart use v=Q_1 and p_old=P_1/a^3.
For the gamma chart the actual canonical swap is

v=P_b/(2*a^3*q), p_old=-2*q*b.

In both charts chi=Q_2 and pchi_old=P_2/a^3. The full
S6.76/S6.77 momentum convention then requires BOTH

p_current=p_old+3*ell*chi,
pchi_current=pchi_old+3*ell*v.

The time-dependent gamma generator was already included
in the scalar Hamiltonian used to derive alpha,beta,V.
All the present phase components fit 2*10^12*U0^2.
There is no Hubble or Theta inversion in this last map.

Passing from globally normalized Fourier columns to
constant-center Fourier coordinates multiplies every
local phase column by a_center^(3/2)<2: the momentum
measure changes by a_center^3, and the normalized
annihilation coefficient changes by its square root.
Both tensor polarization tensors and their duals have
component factors at most two. Reserving a factor four
here gives 8*10^12*U0^2 for all scalar/tensor convention
changes. The common raw seed is much larger:

S=10^16*U0^2=10^42.

Fourier derivatives of any subset of at most four legs
are bounded by K=8U0. Every inverted nonempty proper
transfer has |k_physical|^-1<=D=2/L0. These are the inputs
to the full raw H3/H4 majorant, not a pointwise sampled
replacement for the mode normalization.
