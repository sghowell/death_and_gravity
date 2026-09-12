# Homogeneous Hilbert lift and complete response norm

Work on I=[-1/2,1/2] with smooth tracefree sources compactly supported in its interior, having the common zero initial neighborhood, and Schwartz in space. The unitary spatial Fourier convention is unchanged. Write all tensor norms as Frobenius unless explicitly marked operator norm. Define

    M[D]^2=||D||L2^2+||grad D||L2^2,
    Z136[Gamma]^2=sum_r0..13||(1-Delta)^3 partial_t^r Gamma||L2^2.

S194 bounds the complete homogeneous coordinate-density first current derivative by1e95 for max_j0..12 sup_t||Gamma^(j)||op<=1 and detector ||D||op<=1. A Frobenius-unit detector has operator norm at most one, so finite-dimensional tensor duality bounds the tracefree output's Frobenius norm by the same constant. There is no factor from the five-component tensor basis.

For each real source Fourier component and j<=12, the common zero germ gives

    Gamma^(j)(t)=integral_left^t Gamma^(j+1)(s)ds,
    sup_t||Gamma^(j)||F <=||Gamma^(j+1)||L2_time,F.

The interval has length one. Real linearity of the tangent permits complexification. Apply the bound separately to the real and imaginary source tensors, then add the squared output norms. Their maxima are bounded by the sum of the real and imaginary squared derivative norms through13. Thus no extra sqrt2 factor is necessary. Integrating the output in time and applying spatial Plancherel yields

    |H0[D,Gamma]| <=1e95||D||L2 Z136[Gamma].

This relies on the explicit homogeneous kernel, not an unjustified evaluation atP0 of an arbitrary weakly bounded Fourier operator. The pointwise estimates underlying S208 specialize atP0 before Fourier Cauchy-Schwarz and therefore bound the multiplier Known(0) on the same spatial source space.

Let N61 and X46 be the unchanged S208 norms, Y^2=N61^2+X46^2, and Z24 the S212 two-time/four-spatial norm. For v=|P|^2>=0, the exact nonnegative-polynomial identities for (1+v)^6-1, (1+v)^6-v and (1+v)^6-(1+v)^4 give

    N61<=sqrt2 Z136, X46<=Z136,
    Y<=sqrt3 Z136<2Z136, Z24<=Z136.

S208 bounds each Known term by5e48 M Y; S212 bounds Ffinite by1e5||D||L2 Z24. Hence the unrounded coefficient is at most

    1e95+4*(5e48)+1e5 <2e95.

This proves the full assembled tracefree spatial Gaussian response bound

    |Jren[D,Gamma]| <2e95 M[D]Z136[Gamma]

for nonzero norm product, with equality zero when either source or detector is zero. Extension to the corresponding completion follows by boundedness.

The target is the dual of L2_time H1_space. There is a loss of thirteen time and six spatial derivatives; this is not a same-space contraction. Dividing coordinate density by a^3>=1 cannot increase M of the detector because a depends only on time. Multiplication by both external canonical metric factors gives4/kappa and the display8e-705. None of this supplies a reduced scalar/mixed inverse, stability, or a controlled finite-coupling solution.
