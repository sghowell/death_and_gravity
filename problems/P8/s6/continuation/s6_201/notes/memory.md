# Full time-ordered memory correction

The unit-W8 full stress-pair coefficient satisfies the existing
Cref=1e9 bound Cref sqrt(nu mu). This is conservative: the direct
unit-readout count is 2*4000^2, without either alpha enlargement.

The exact occupation increment therefore bounds each detector/source
product by

    2 Cref^2 B^2 nu mu(nu^-12+mu^-12)
      ||Dhat(t,P)||F ||Gammahat(s,P)||F.

Retain theta(t-s) and all nine physical polarization pairs.
For P=k+l set L=1+|P|/(Amax m), Amax=25/16. The elementary massive
frequency comparison gives mu<=L nu and nu<=L mu. Hence

    integral nu mu(nu^-12+mu^-12) <=2L J10,
    J10=5 Amax^3/(512pi m^7).

The exact full pair factor is thus 36 Cref^2 B^2 L J10.
The unit time triangle is bounded by time L2 Cauchy without
differentiating or altering theta. Since L<2(1+|P|^2), spatial
Plancherel and weighted Cauchy give

    |R_alpha,memory| <=72 Cref^2 B^2 J10 M[D]M[Gamma]
                      <1e12 M[D]M[Gamma]

for nonzero input norms. This is an absolute continuum bound.
No initial mixing derivative, spectral reset, discarded quadratic
occupation term or restriction of external momentum is used.
