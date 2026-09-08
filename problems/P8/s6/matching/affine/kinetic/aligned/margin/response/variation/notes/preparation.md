# The actual response requires a changed initial covariance

The issue is not an arbitrary hypothetical metric variation.
For the S6.57 zero-independent-data response, the initial
lapse is n0=-F_n(u0)/(2J_e(u0)), u0=-1/2.

At that slice the exact local force coefficients in units
R^(4-2j)/(64pi²L²), j=0,1,2, are

    32651/6750, 13805696/253125, -183496192/6328125.

Their sum at R=1000 exceeds 4R^4. Using pi<22/7 and
(22/7)²<10 gives a local lower bound R^4/(160L²).
The full subtracted-state error, including the S6.55 state
change, is at most

    6429601/(72000000000000000000000000).

At L=10^12 the resulting F_n lower bound is above
5*10^-15. Also 0<J_e(u0)<1 for the entire allowed
epsilon range. Therefore n0<-1/(4*10^14), while its
absolute value is below S6.57's 94*10^-14 upper bound.
In particular 0<N0=1+n0<1 and the original clock tube
is retained. The negative curvature-order-four term
and the nonzero state error are both included.

## Physical longitudinal covariance

At the initial slice v=0, so the exact physical scale ratio
is exp(omega). For the high-k physical canonical oscillator,
the relevant positive-frequency impedance is g² Omega.
For transverse modes it tends to k_com, independent of
a and N. This is a conformal leading-order control.

For longitudinal modes it instead tends to

    a² m² k_com sqrt(a_m*b_m).

Relative to the old initial geometry, its ratio is
r_L=exp(2omega)sqrt(a_m*b_m). Since
exp(2omega)=1/(2p_affine),

    r_L²=a_m*b_m/(4p_affine²).

N0<1 implies p_affine>1/2. The original tube gives
p_affine<21/40. From the unchanged exact mass functions,

    r_L²-1=-(2p-1)Q(p)/
              [324p²(8p+5)(2p³-1)],
    Q=2592p^5+3708p^4+886p³-1350p²-2013p-605.

On [1/2,21/40], bound positive powers at 21/40 and
negative powers at 1/2. This proves Q<0. The mass
numerators/denominators keep a_m,b_m positive and
2p³-1<0. Consequently 0<r_L²<1 at the actual
initial lapse. The on-clock logarithmic slope is the
independently checked 113/(81h), not zero.

For old and new leading impedances Z and rZ, copying
the old physical coordinate/momentum data produces

    A_infinity=(r+1)/(2sqrt(r)),
    B_infinity=(r-1)/(2sqrt(r)),
    |A_infinity|²-|B_infinity|²=1.

The copied data preserve CCR but have a nonvanishing
negative-frequency principal component when r!=1.
Their vacuum-subtracted oscillator excitation energy
has a positive k*|B_infinity|² leading term. The radial
integral then diverges like integral k³ dk. Smallness
of n0 does not make this UV tail integrable.

Thus blindly copying the old longitudinal physical
Cauchy covariance fails even the necessary leading
adiabatic preparation condition for the new operator.
This does not exclude a compatible all-order state family.

## Different canonical prescriptions are not interchangeable

Freezing the normalized oscillator data is a different
prescription from freezing the physical canonical pair.
Even the transverse oscillator leading frequency then
changes by the factor r_T=N*exp(-omega), with

    r_T^4-1=(N²-1)[(h-1)N²+h]/h.

For positive N and h>=1, r_T=1 only at N=1. Hence that
different frozen-data prescription already fails at
leading transverse order. The physical transverse
conformal control does not establish all its subleading
Hadamard data.

Off the exact clock the temporal and spatial mass
functions differ, and the longitudinal principal
operator is not ordinary Proca in the same physical
metric. S6.55's ordinary-Proca propagation theorem
is not silently applied to that different operator.

The required next step is a compatible initial covariance
family and its finite, renormalized causal response.
The existing fixed-source result remains valid within
its stated boundary. No change of state or new input
from the user is being presumed here.
