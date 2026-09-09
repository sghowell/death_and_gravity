# Physical boosts, actual profile Gram and two nonzero widths

Let partial_plus=(partial_t+partial_z)/2 and
partial_minus=(partial_t-partial_z)/2. A future unit timelike and
orthogonal unit spacelike pair is

    U=b partial_plus+b^-1 partial_minus,
    E=b partial_plus-b^-1 partial_minus, b>0.

Its null derivative U+E=b D. The flat physical stress contraction
in (1) must therefore be multiplied by b^-2 to bound T_DD.
The Lorentz determinant is one. The state term is invariant because
b^-2[(U+E)f]^2=(Df)^2.

For f=sqrt(2)h(x_plus)j(x_minus) with each one-dimensional norm one,
the plane measure is dx_plus dx_minus/2, so ||f||=1. Write
H1=||h'||^2,H2=||h''||^2 and J1,J2 similarly. Compact integration
by parts gives <h,h'>=<h',h''>=0 and <h,h''>=-H1. Consequently,
in the ordered basis (f_++,f_+-,f_--), the derivative Gram is

    [ H2      0      H1 J1 ]
    [  0    H1 J1      0   ]
    [ H1 J1   0        J2  ].                             (10)

Using U^2 coefficients (b^2,2,b^-2), and UE coefficients
(b^2,0,-b^-2), the three derivative pairings are

    ||U^2 f||^2=b^4 H2+6H1 J1+b^-4 J2,
    ||UE f||^2=b^4 H2-2H1 J1+b^-4 J2,
    <U^2 f,UE f>=b^4 H2-b^-4 J2.

Thus the quantum cost, before hbar/(8pi^2), is

    (2/3)(1+4xi)b^2 H2+4(1-2xi)b^-2 H1 J1
                       +(2/3)b^-6 J2.                   (11)

The state cost is 8xi Phi_*^2 H1. Rescale h,j by the positive
widths delta_plus,delta_minus and set b^2=r delta_plus/delta_minus.
Every term in (11) then has the common denominator
delta_plus^3 delta_minus. No width is dropped or identified with
a UV cutoff.

For the actual clamped polynomial in the formulation, native
integration and a separate Fraction-only coefficient convolution
give

    <h^(i),h^(j)>_(i,j=0,1,2)
       =[[1,0,-3],[0,3,0],[-3,0,63/2]].                  (12)

The actual product Gram obtained from (12) independently agrees
with (10), including its off-diagonal 9. The profile and its first
derivative have zero outer traces. Its extension by zero belongs
to H0^2(-1,1), although its second derivative jumps at the ends.
Approximate it by C_c^infinity(-1,1) functions in H^2 and normalize.
Their products converge in the required plane derivative norms.
For each fixed Hadamard target the stress and Wick expectations
are smooth and bounded on the compact rectangle, so their weighted
integrals converge as well. The one-sided cap is needed only on
this rectangle; the approximating supports remain inside it.

At xi=1/6, (11) yields Gamma(r)=35r+24/r+21/r^3. Its derivative
vanishes when 35r^4-24r^2-63=0, with the unique positive solution
r^2=(12+9sqrt(29))/35. Gamma''=48/r^3+252/r^5>0 and Gamma diverges
at both ends of r>0. The minimizer lies strictly between 1 and 4/3,
as checked by opposite signs of the stationary polynomial there.
Choosing r=4/3, without approximating that root, gives

    Gamma=14117/192, 74-Gamma=91/192>0,
    state cost=4 Phi_*^2/delta_plus^2.                    (13)

This proves (2). The algebraic minimizer is ONLY a minimizer over
this boost family for this profile. No claim of optimal sampler,
optimal positive-type decomposition or optimal physical lower
bound is made. The two-null-width scaling in (13) explicitly
diverges in the attempted finite-length single-null limit.
