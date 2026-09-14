# Full sixteen-phase nonautonomous cone and physical conversion

Let eplus=(1,1,1,1)/2 and eminus=(1,-1,-1,1)/2 in the four-cycle
block. Its symmetric part is exactly
eplus eplus^T-eminus eminus^T. Both vectors have unit norm; the
positive vector is a right and left eigenvector with eigenvalue 1.
On its orthogonal complement the cycle logarithmic norm is at most
zero. The entire retained slow block has logarithmic norm
d/(2*80000)<=rho_min/8, with rho_min=10^-4. All five remaining
oscillators (H, two TT, two transverse Proca) are retained skew cores
in their energy variables, even when their frequencies are large.

Normalize the equation by P^3/2 without rescaling the time-dependent
coefficients away. The complete matrix remainder and every chart
connection together have norm at most

    error=10^-19+10^-33+10^-46 < delta=rho_min/32.

Split the full sixteen-phase vector into f eplus+z, with z including
the entire orthogonal complement, not merely the other fast roots.
For f>0 and ||z||<=f, the whole equations give

    f'/P^3/2 >= (rho_min-2 delta) f.

On the cone boundary ||z||=f, the upper derivative satisfies

    D+||z||/P^3/2 <= (rho_min/8+2 delta) f.

The factors 2 bound sqrt(2) from the full remainder norm. The gap
between these two boundary coefficients is 3 rho_min/4>0.
Consequently the forward cone is strictly invariant. Start with
f=1,z=0 and integrate the weaker rate rho_min/2. Uniformly over
the ENTIRE finite band, at T=10^-60,

    ||U_growth(T,0)|| >= exp(rho_min P_min^3/2 T/2)
                        =exp(5*10^31).

This is a varying-coefficient full-generator argument; an instantaneous
quartic determinant alone would not prove it. The positive heavy mass
is retained in the complement, not used as a disastrously large raw
remainder estimate.

To return to physical canonical variables, let Perm be the explicit
eight-phase reordering in the matrix note. The full coupled maps are

    Kphysical=W M^-1 Perm R^-1 T0,
    Kphysical^-1=T0^-1 R Perm^T M W^-1.

Both entire matrices are formed and factored before enclosure, so
momentum powers that cancel between transformations actually cancel.
Each operator norm is below 10^150 on the entire band and coefficient
box. The exact product is the identity. This includes the physical
central transformation AND the perfect-square shear, not just the
final energy weights. Every TT/transverse-Proca energy chart and
inverse independently fits the same bound. Their direct sum has
operator norm equal to the largest block norm, so the full sixteen
phase maps and inverses also obey 10^150.

It follows, by choosing the initial physical vector Kphysical(0)^-1
eplus and using the terminal forward-map bound, that

    ||U_physical(T,0)|| >= 10^-300 exp(5*10^31).

The common sqrt(kappa) canonical normalization cancels in this
propagator comparison; no physical field is rescaled differently
between initial and final times.

All interval calculations enclose every real P in the band; the
exact-rational input validator is a certificate interface, not a
restriction of the continuum theorem to rational momenta. Smooth
Fourier packets in a small directional sector of the band, paired
with their real conjugates, have the same L2 lower bound by Parseval.
The initial maps and polarizations can be chosen smoothly there.
No limit P->infinity is used.

The actual metric relation gives a_phys=a_hat/R^1/4. The full
coefficient box encloses R^1/4/a_hat between 1/2 and 2, hence
physical wavenumbers lie in [10^64/2,4*10^64]. Their square is
strictly below n. Also rho P^3/2<3*10^94 throughout the band;
proper time divides by N, so even 6*10^94 has square below n.
The retained heavy oscillator is not claimed to have frequency
below its own mass. These sub-heavy comparisons do NOT establish
the actual Wilsonian cutoff or bound errors from a different
matched EFT or UV completion.
