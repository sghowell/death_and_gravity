# Global physical matter observable, not a fictitious fluid

Use the actual intrinsic physical density and pressure
Hessians derived in S6.106, including inhomogeneous lapse,
scalar matter charge and spatial-gradient fluctuations.
Transform them by B and T into W coordinates.
Multiply each kernel by t^4 and apply the same exact
pointwise envelopes with one extra half-order for each
momentum index. Every resulting entry is uniformly bounded.

For each physical observable the code sums HALF of all
entry bounds, including off-diagonal entries, and multiplies
by C_cov. The homogeneous mean contribution is counted once:

    delta rho_mean=delta p_mean=-3ell^2[xi+n/(2h)].

Since ell=1/(10t^6), its magnitude is at most
(3/100)C_B eta_S/t^12. Combining intrinsic and mean terms,

    |delta rho_m|, |delta p_m|
       <10^64 eta_S/t^4

in physical units M^2/tau^2. This bound uses the actual
mean lapse after trace restoration. No second mean lapse
is inserted inside the intrinsic quadratic term.

The observable is the coordinate spatial average of the
clock-normal density/pressure jet. It is not an independently
conserved homogeneous fluid, total gravitational tadpole,
or perturbed-volume integrated energy. The parent center
comparison to the independently constrained physical density
remains unchanged by the analysis-variable compensation.

The bound is absolute, NOT uniform fractional smallness
relative to the classical matter density ell^2/2~t^-12.
An upper bound decaying as t^-4 cannot establish such
fractional smallness. Nor does that upper bound by itself
prove that the actual fractional response diverges.
No QEI, positivity, full Ward or higher-order conclusion
is inferred from a sign or the decay rate of this jet.
