# Continuous all-momentum adjoint wave estimate

All L2 norms use dt dx and the spatial TT Frobenius norm.
For a real vector of components solving Lv=q, define

    E(t)=1/2 integral[a^3 |v_t|^2+a |grad v|^2+a^3 |v|^2]dx.

Spatial integration by parts in the actual wave equation
gives exactly

    E'=integral a^3 v_t(q+v)dx
       +(H/2)integral[-3a^3 |v_t|^2+a |grad v|^2+3a^3 |v|^2]dx.

The cross term is bounded by E and the H terms by
3|H|E. Cauchy-Schwarz controls the forcing. Therefore

    |E'| <= (1+3|H|)E+sqrt(2E)||a^(3/2)q||L2(dx).

On this actual slab, 1<=a<=25/16, |H|<=8/5 and
|H'|<=4. For y=sqrt(2E), regularizing E by a positive
epsilon and then taking its limit justifies the estimate
also at zero energy:

    |y'| <= 3y+||a^(3/2)q||L2(dx).

Integrating backward from zero terminal data, using length1,
exp(3)<27 and (25/16)^(3/2)=125/64<2, gives

    sup_t y(t)<54 ||q||L2(dt dx).

Consequently the spacetime norms of v, v_t and grad v
are each below54 ||q||. The added v^2 energy term
controls zero momentum: no mass, Fourier band or inverse
momentum factor has been introduced.

Define the detector norm

    D(q)^2=sum_(|alpha|<=2)||partial_x^alpha q||L2^2
           +||partial_t q||L2^2,

where alpha ranges over all spatial multiindices.
The exact time-dependent coefficients commute with Delta.
Since ||Delta q||<=sqrt(3) D(q)<2D(q), apply the same
energy estimate to Delta v to get
||Delta v||,||Delta v_t||<108D(q).

From the actual wave equation,

    v_tt=q-3H v_t+a^-2 Delta v,

the spacetime norm is below
[1+(24/5)54+108]D=368.2D<400D. Differentiating it once,

    v_ttt=q_t-3H v_tt-3H'v_t
           +a^-2 Delta v_t-2H a^-2 Delta v,

gives a bound by
[1+(24/5)400+12*54+108+(16/5)108]D
=3022.6D<4000D. These are full spatial estimates.

Thus the complete S6.186 stress-test norm satisfies

    N(v)^2 <= [3*54^2+400^2+4000^2]D(q)^2
            =16168748 D(q)^2 <5000^2 D(q)^2.

The energy estimate holds for all compact q and its
advanced solution on the closed slab. Applying the
stress-distribution bound to an extension across its
initial endpoint requires the two compatibilities in
boundary.md. The estimate alone does not supply them.
Complex Fourier components follow by applying the same
real estimate to real and imaginary parts, or directly
using Hermitian inner products.
