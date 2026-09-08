# Off-clock principal-speed check for the retained vector block

S6.58 derives the full temporal elimination for homogeneous physical
histories with N,a,a_mass,b_mass>0. Its bare canonical frequencies are

    Omega_T^2=N^2(q+m^2 b_mass),
    Omega_L^2=N^2(b_mass*q/a_mass+m^2 b_mass), q=k_com^2/a^2.

Canonical normalization subtracts d'+d^2. With H=a'/a,

    d_T=(H-N'/N)/2,
    d_L=[H+a_mass'/a_mass-N'/N
           +(2Hq-m^2 a_mass')/(q+m^2 a_mass)]/2.

The code differentiates these full rates, including q'=-2Hq,
and verifies that (d'+d^2)/q tends to zero. Normalizing the
large-q coefficient by the physical metric's N^2 gives

    c_T^2=1, c_L^2=b_mass/a_mass.

These are principal speeds of the retained two-derivative vector
block on a prescribed metric. They are not full coupled, quantum
or UV signal-front verdicts.

For the actual unchanged mass functions,

    c_L^2-1=(2p_affine-1) F(p_affine),
    F=-3p(24p^3+12p^2+6p-5)/[(8p+5)(22p^3+3p-11)].

On 9/20<=p<=21/40 the cubic in the numerator is positive,
the temporal-mass numerator is negative, and both masses
are positive. Endpoint monomial bounds prove all signs on
the entire interval. Thus F>0. This interval contains the
original open clock tube p^2 in [9/40,11/40].

The original point chart has

    4p^2-1=(N^-2-1)/h, h=(1+u^2)^3>0.

Hence this block is strictly superluminal for N<1, luminal
at N=1 and subluminal for N>1. At the clock its lapse slope
is -8/(81h), obtained from dc_L^2/dp=16/81 and
dp/dN=-1/(2h).

S6.58 proves that the particular zero-independent-data S6.57
fixed-source response has N_initial<1. The retained vector
block on that approximate metric therefore fails this cone
test at its initial slice. At fixed m0*tau=1000, both sides
of the source estimate scale with (M*tau)^-2, so increasing
M*tau decreases the excess without changing its sign.

This neither excludes another initial-data choice nor establishes
the complete coupled quantum cone: mixing, quantum corrections
and omitted operators remain to be bounded. It does show why
the old approximate metric should not simply be promoted to
a certified causal quantum bounce.
