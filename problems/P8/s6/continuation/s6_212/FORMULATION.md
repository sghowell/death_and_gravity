# Formulation: original dimensional spatial UV finite difference

## Fixed problem and new calculation

Retain the original CD slab |t|<=1/2, a=(1+t^2)^2, H=a'/a, m=1000, kappa=10^800, ell=log(m^2/mu^2)=0, original physical metric, prepared state, unit-W8 comparison and fixed subtraction prescription.

Continue only the singular mode-expression calculation to d=3-2epsilon spatial dimensions. The already convergent actual-state difference stays in physical dimension3. No noninteger-dimensional physical Hilbert space is asserted.

For tracefree detector/source tensors and P=p e, use T=tr(D Gamma), V=(De).(Gamma e), W=(e.De)(e.Gamma.e). Normalize the sphere averages before restoring the full Fourier factor. Reconstruct invariant T,V,W coefficients before differentiating d: the normalized scalar test tensor otherwise changes the source.

## Exact finite coefficient

Let F2 and F4 be the physical complete spatial UV differences from S211, including every original endpoint and source-time slot. Let f4(d) be the normalized-sphere complete logarithmic difference for the fixed source invariants. Let Hpole(d) be the full spatial second variation of

    m^2 R_old-R_old^2/4+29Ricci^2/30-2Riemann^2/15,

with these fixed four-dimensional scalar weights but continued metric contractions and volume. The local spatial UV finite difference is

    Ffinite = -m^2 F2 Gamma/2
              +(1-log2-ell/2)F4 Gamma
              -(partial_d f4)(3)/(2pi^2)
              +(partial_d Hpole)(3)/(32pi^2).

The source-value fourth-grade term includes j1/r0/d3. Its physical value is zero but its dimension derivative is nonzero. Neither this contribution nor the continued Euler and volume terms can be discarded before taking the finite limit.

The radial normalization is the original MSbar one. The split begins at the fixed comoving radiusm, not a*m, so the lower finite polynomial and volume logarithm remain. No finite counterterm or renormalization scale is selected anew.

The first-time coefficient of Ffinite is exactly the proper-time derivative of its second-time coefficient. An explicit log-polynomial bound gives

    |Ffinite[D,Gamma]| < 1e5 ||D||L2 Z24[Gamma]

for nonzero finite-norm smears, with the zero case understood, where Z24 has up to two source-time and four spatial Sobolev derivatives. Both canonical metric factors give4e-795 for this coefficient alone.

## Boundary

The full dimension-limit interchange and assembly with the fixed homogeneous anchor and known actual remainders still require proof. No complete spatial response/inverse, mixed scalar normalization, finite-coupling background/stability, remaining parent loop/cutoff matching, vacuum/finite-gravity gate or original P8 closure is claimed.
