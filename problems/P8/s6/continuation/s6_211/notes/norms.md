# Explicit coefficient bounds, not cutoff-uniform response bounds

Use the unchanged Fourier measure dt d^3P/(2pi)^3 and

    Z24[Gamma]^2=sum_(r=0)^2 integral
       (1+|P|^2)^4 |partial_t^r Gamma_hat|_F^2.

For unit Frobenius tensors, |T|,|V|,|W|<=1; hence |T-2V|<=3 and |T-V|<=2. On the original slab,1<=a<=25/16, |a'|<=5/2 and aQ=4+28t^2<=11. Use pi^2>9.

The absolute coefficient row for DeltaF2 is at most97 p^2/(3360*9), and97/30240<1/200. The five DeltaF4 row bounds, weighted by(1+p^2)^2, are respectively

    mass: 390625/48,
    curvature: 11/576,
    spatial: 17/5760,
    source first: 13/1728,
    source second: 65/13824.

Their exact sum562502369/69120 is below10000. Plancherel followed by spacetime Cauchy-Schwarz gives

    |DeltaF2[D,Gamma]| < (1/200)||D||L2 Z24[Gamma],
    |DeltaF4[D,Gamma]| < 1e4 ||D||L2 Z24[Gamma]

for nonzero finite-norm smears, with the corresponding non-strict zero bound for zero inputs. Polynomial symbols ensure the P0 extension is regular despite the direction notation.

Both canonical metric insertions multiply these bounds by4/kappa; the logarithmic coefficient display is4e-796. The K^2 and logK factors still grow. This is not a bound for a fully matched response, not a canonical scalar or full mixed norm, and not a controlled inverse, nonlinear solution or stability estimate.
