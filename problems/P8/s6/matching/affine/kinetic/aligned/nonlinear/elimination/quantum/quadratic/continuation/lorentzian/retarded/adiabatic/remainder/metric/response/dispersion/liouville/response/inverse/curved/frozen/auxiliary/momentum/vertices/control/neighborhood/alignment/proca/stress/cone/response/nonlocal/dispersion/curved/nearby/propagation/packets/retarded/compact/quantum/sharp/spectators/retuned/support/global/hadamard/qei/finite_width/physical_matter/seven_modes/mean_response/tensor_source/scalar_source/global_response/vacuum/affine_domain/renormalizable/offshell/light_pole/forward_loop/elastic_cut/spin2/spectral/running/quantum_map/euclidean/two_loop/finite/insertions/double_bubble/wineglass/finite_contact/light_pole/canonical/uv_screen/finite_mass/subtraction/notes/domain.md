# A complex domain for the first-cut spectral function

Keep Q=16pi^2, a=g_s^2, Y=y^2 and m=mF.
For |s-2|<=5, both |s| and |s-4| are at most seven.
The on-shell parameterization E^2=s/4, p^2=s/4-1
therefore has |E|,|p|<=sqrt(7)/2<3/2.
Use the same real gauge directions and transverse
polarizations as S6.129.

Each Euclidean scalar vector has component norm below
three. Each gauge vector has norm at most
|E|(1+sqrt(3))<9/2. Every one of the six cyclic routes
contains at most three cumulative external vectors,
so every propagator shift has norm below 27/2<18.
For m>=36, the strict Neumann ratio is below 18/m<=1/2.

The complete regulated Ward sum and the leading
coefficient are unchanged. Applying the already proved
all-degree radial/composition majorant with routing
constant 18 instead of 12 gives

    |R_box| < 96*6*18^3 aY/(Qm^3)
            =3359232 aY/(Qm^3)
            <4*10^6 aY/(Qm^3).

This is the whole one-loop momentum tail, not a finite
sample of expansion terms or an all-loop remainder.

## Descent through kinematic square roots

The box sum is holomorphic in the external complex
components on this uniform domain. Interchanging the
identical scalar legs sends p to -p and leaves the
complete transition unchanged. A rank-two Lorentz
tensor built from external momenta has even total
momentum degree; equivalently the simultaneous
(E,p)->(-E,-p) inversion is even when the two
polarization vectors are kept fixed. There is no
Lorentz-invariant odd-rank coefficient. Consequently
the transition is separately even in E and p.

The same holds for the product of the two continued
transition factors. Its convergent power series
descends to E^2 and p^2. This removes the artificial
square-root monodromies at s=0 and s=4. Compact real
angular integration preserves holomorphy and the
uniform bound on the whole stated s disc.
It does not conjugate an unphysical complex momentum.
No new heavy singularity or uncontrolled complex
loop-contour translation is introduced.

The analytic off-shell Ward identities also give
a soft factor on each gauge leg. For example, expand
the tensor in independent k near zero, holding l
and the independent scalar momentum fixed. The
coefficient of k in k_mu T_mu_nu=0 forces the
k-independent tensor to vanish. Repeat with l.
Thus B=O(k l)=O(E^2) on the soft on-shell slice, and
rho=O(s^2) at the massless cut endpoint. This is
a statement about this massive one-loop transition,
not a general all-loop infrared theorem.

## Uniform continued cut error

Let z=aY/(Qm^2), epsilon=4*10^6/m and
eta=5epsilon/2+epsilon^2/8. On the complex s disc
the leading pair has absolute magnitude at most
(4/3)z|s|<=28z/3<10z. Each correction is bounded
by epsilon z. Summing all four pair products gives
(80epsilon+4epsilon^2)z^2. The unchanged phase,
optical, identical-particle and color factors yield

    rho_0(s)=dA z^2 s^2/(9pi)=pi K s^2,
    |delta rho(s)|<=M_rho=dA z^2 eta/pi=9pi K eta,
    K=dA z^2/(9pi^2)=dA C^2/pi^2, dA=8.

The continued rho and delta rho are holomorphic on
|s-2|<=5. Positivity of a physical massive-scalar
cross section below threshold is not assumed.
