# Full angular dimension jets at a fixed physical source

All105 complete angular channel/source/inverse-radius coefficients are evaluated with the actual dimensional modes and geometry. Atd3 they reproduce every full S211 coefficient, not only the nonzero spatial differences.

For each coefficient, reconstruct the invariant row before differentiating:

    A=f_tensor,
    B=2(f_vector-f_tensor),
    C=d(f_scalar-f_tensor)/(d-1)-B.

The scalar probe has V=W=(d-1)/d. Differentiating its raw channel coefficient atd3 would therefore include an unintended source variation. The exact difference from the fixed-source derivative is(B+C)(3)/9. This is checked for every slot and is nonzero in relevant fourth-grade coefficients.

The complete logarithmic source coefficient is the sum over all original endpoints j and inverse-radius degree4-j. The j1/source-value/degree3 spatial difference vanishes atd3, but its dimension derivative does not. Explicit invariant derivatives of this endpoint are retained in the certificate. Higher source slots and finite physical odd endpoints are not deleted.

Tables use normalized sphere averages. The full physical Fourier factor is

    M(d)=2pi^(d/2)/[Gamma(d/2)(2pi)^d],
    M(3)=1/(2pi^2),
    partial_d log M(3)=EulerGamma/2-1-log(pi)/2.

In d=3-2epsilon, both M(d) and the coefficient derivative contribute finite terms. The analytic formula is not obtained by freezing the number of transverse polarizations or the three-dimensional magnetic dual.

The one-leg contact remains local in the pointwise metric matrices for general d. For external P,-P its convolution has total0 before the homogeneous covariance trace, hence its spatial difference vanishes before internal integration. This retains the generally nonzero homogeneous contact and does not introduce a new state/history.
