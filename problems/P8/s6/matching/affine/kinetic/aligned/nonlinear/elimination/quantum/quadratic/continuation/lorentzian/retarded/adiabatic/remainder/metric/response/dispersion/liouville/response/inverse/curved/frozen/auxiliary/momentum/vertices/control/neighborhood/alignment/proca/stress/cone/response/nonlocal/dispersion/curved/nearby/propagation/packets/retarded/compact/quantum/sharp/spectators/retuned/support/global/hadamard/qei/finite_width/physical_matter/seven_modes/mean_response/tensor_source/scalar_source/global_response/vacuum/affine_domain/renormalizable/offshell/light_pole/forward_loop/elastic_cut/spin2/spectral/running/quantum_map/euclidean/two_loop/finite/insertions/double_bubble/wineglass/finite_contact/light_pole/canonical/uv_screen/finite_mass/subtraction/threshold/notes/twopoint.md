# Independent two-point trace and complex pole remainder

With Hermitian Euclidean gamma matrices, tr I=4, the
numerator is

    tr[(mF-i gamma.q)(mF-i gamma.(q+p))]
      =4[mF^2-q.(q+p)]
      =2(4mF^2+pE^2)-2(Dq+Dq+p).

The Clifford trace reduction is valid for the regulated
contractions with tr I=4; the explicit matrix test checks
its finite algebraic numerator. We do not replace an
unsubtracted D-dimensional integral with a four-dimensional
norm estimate. Shift/reduction identities are used in the
same translation-invariant dimensional regulator.

Expansion of minus the fermion log determinant to second
order gives 2NY[(4mF^2+pE^2)B-2T]. Continuing s=-pE^2,
at nu=mF the finite tadpole is -mF^2/Q and the bubble
is J(s)/Q, giving the f(s) in the formulation.
The entire ultraviolet affine reference is
2NY(6mF^2-s)Ibar/Q. Subtracting at s=1 through first
derivative cancels it exactly, together with the
constant tadpole. The finite mass reference remains
explicit; it is not discarded by the remainder operation.

For |s-1|<=2, |s|<=3 and Re(mF^2-As)>=mF^2-3/4>0.
There is no parameter pinch; the logarithm and its
derivatives are holomorphic. The second derivative of
(4mF^2-s)[-log(1-As/mF^2)] is

    [mF^2(4A^2-2A)+s A^2]/(mF^2-As)^2.

Bound its absolute numerator by mF^2(4A^2+2A)+3A^2.
The moments integral A=1/6 and integral A^2=1/30 give
(7/15)mF^2+1/10. Taylor's integral remainder along the
straight segment from 1 to s supplies a factor 1/2,
canceling the leading factor 2 in f. This proves
|f_R(s)|<=B|s-1|^2 uniformly, including complex s.

On 0<=s<=1 the second derivative is strictly negative:
the numerator is at most -A(mF^2-sA)<0 away from
the endpoints. Hence f'(1)<f'(0)=4NY/(3Q); the series
lower endpoint proves positivity at the named mass.

For s>4mF^2, the physical boundary of J has imaginary
part pi sqrt(1-4mF^2/s). Consequently
Im f=-NY s(1-4mF^2/s)^(3/2)/(8pi).
Our continued Euclidean inverse is m_ref^2-s+f;
its negative is the usual Minkowski inverse with the
positive absorptive width. Do not import a self-energy
sign convention from a different inverse propagator.
