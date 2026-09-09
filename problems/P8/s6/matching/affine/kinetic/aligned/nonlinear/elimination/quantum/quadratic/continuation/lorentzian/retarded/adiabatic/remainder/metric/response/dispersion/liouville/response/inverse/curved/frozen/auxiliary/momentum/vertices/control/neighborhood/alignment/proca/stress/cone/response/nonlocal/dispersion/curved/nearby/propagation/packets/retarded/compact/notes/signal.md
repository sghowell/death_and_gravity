# Positive compact response with the entire L2 remainder retained

Let a,S_min,S_max and A_min be as in notes/probes.md.
For |t-t0|<=a/2 and any source point in supp J, take sphere
directions n with n1>0 and |n_perp|<=a/(4S_max).

The source shifts the clock radius from S_c(t,s0) by at most
4rho and its spatial center by at most rho. On
x=y+S_c(t,s)n the detector normal displacement is therefore
at most 5rho=epsilon/4. Its transverse displacement is at
most rho+a/4<a/2. Its longitudinal displacement is at most

    5rho+a^2/(16S_max)<a/2,

using 1-sqrt(1-z)<=z on [0,1]. All four detector bumps
are at least one half on this cap, so f>=1/16.

The cap solid angle satisfies

    2pi(1-sqrt(1-b^2))>=pi b^2,
    b=a/(4S_max).

The spatial integral of the clock shell is
A_c S_c/(4pi) times the angular integral. Integrating the
time core of length a and the normalized source yields
the exact strictly positive lower bound

    D=A_min S_min a^3/(1024 S_max^2).

No source normalizer, sphere radius or time width is discarded.

## Detector norm and exact remainder size

For fixed t the transverse support is a disk of radius a.
The radial band has thickness 2epsilon. It lies on the
positive x1 branch, away from the origin. With
a<S_min/4 and epsilon<S_min/4, |x_perp|/r<=1/3, hence
dx1/dr<2. The support volume is at most

    4pi a^2 epsilon <16 a^2 epsilon.

Since |f|<=1 and its time support has length 2a,

    ||f||_(L1_time L2_space) <=8a^2 sqrt(epsilon).

S6.91 gives a uniform spatial L2 remainder. Write
K0=4*10^31, C=10^44 and let B be its full low-frequency
original-generator bound. Define

    M_low=16 exp(T B)+8*10^6 T,
    R_bar=C+K0^2 M_low.

This bounds the complete remainder norm: the high piece
has norm at most C/(sqrt(2) pi sqrt(K0))<C; the low piece
has norm at most K0^(3/2) M_low/(sqrt(6) pi)<K0^2 M_low.
The triangle inequality combines both pieces. This
overbound is finite and greater than C, but extremely
conservative. No decimal evaluation of exp(T B) is used.

Choose the exact finite widths

    E0=(D/(32a^2))^2,
    epsilon=E0/R_bar^2, rho=epsilon/20.

The rational bound epsilon<E0/C^2<a/100 verifies every
geometric restriction above. Since the source has spacetime
L1 norm one, translation invariance and Cauchy-Schwarz give

    |remainder pairing|<=R_bar ||f||_(L1 L2)
                        <=8R_bar a^2 sqrt(epsilon)=D/4.

The matter shell is absent and the clock shell contributes
at least D. The actual full real compact pairing is therefore

    C_classical >=3D/4 >0.

This is a quantitative lower bound for specified mathematical
probes. It is not a detector noise/resolution bound, nonlinear
solution or frequency range admitted by an interacting EFT.
The tiny width is retained explicitly rather than hidden in
an unspecified localization limit.
