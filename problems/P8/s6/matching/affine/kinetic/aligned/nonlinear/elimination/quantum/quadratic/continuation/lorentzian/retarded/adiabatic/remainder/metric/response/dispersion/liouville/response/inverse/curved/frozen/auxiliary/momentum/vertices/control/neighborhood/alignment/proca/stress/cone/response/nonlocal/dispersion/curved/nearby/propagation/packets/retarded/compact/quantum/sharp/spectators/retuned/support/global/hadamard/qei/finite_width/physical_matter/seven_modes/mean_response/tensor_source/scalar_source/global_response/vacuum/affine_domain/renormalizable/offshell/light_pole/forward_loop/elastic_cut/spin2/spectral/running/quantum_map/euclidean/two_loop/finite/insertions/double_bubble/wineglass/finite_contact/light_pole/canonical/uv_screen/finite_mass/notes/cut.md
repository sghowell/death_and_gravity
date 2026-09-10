# Finite-mass first two-gauge discontinuity

Write z=aY/(Qm^2)>0 and epsilon=10^6/m.
The leading local operator has C=-z/3 and vertex 8 C T.
The inherited independent polarization/Lorentz contraction
gives two parallel amplitudes of magnitude 4zs/3 and
two orthogonal amplitudes zero. For 1<=s<=3 every leading
pair is bounded by 4z, and every finite-mass correction
by epsilon z, uniformly over the cut angle.

For both holomorphically continued transition factors,
the elementary product estimate is

    |(B0+R)(B0_tilde+R_tilde)-B0 B0_tilde|
      <= |B0| |R_tilde| + |B0_tilde| |R| + |R| |R_tilde|.

There is no assumption that the two complex factors
are conjugates below the scalar threshold. Summing all
four transverse pairs gives at most
(32 epsilon+4 epsilon^2) z^2. The color sum is dA=8.
Integrated two-massless-particle phase space is 1/(8pi),
the identical-particle factor is 1/2 and the optical
conversion 2 Im A supplies another 1/2. Therefore

    rho_s^0 = dA z^2 s^2/(9pi),
    |rho_s-rho_s^0| <= dA z^2/pi
                               (epsilon+epsilon^2/8).

Here rho denotes the continued single-channel cut divided
by 2i. It is not asserted to be a positive cross section
in 0<s<4. At forward transfer u=4-s, the upper and lower
physical s boundaries give opposite u boundary signs.
The total boundary jump is 2i(rho_s-rho_u).

At s=3/2 and u=5/2 the leading difference equals
-4 dA z^2/(9pi). The sum of the two remainder bounds
gives a relative error at most

    (9/2)(epsilon+epsilon^2/8) < 5 epsilon,

provided 0<epsilon<=1/2. The strict margin is
epsilon(8-9epsilon)/16>0. For the unchanged m=10^200,
epsilon=10^-194, so the exact finite-mass first cut
cannot cancel this difference at its perturbative order.

The general rational interface keeps the leading
normalized difference D=s^2-(4-s)^2=8(s-2) and error
18(epsilon+epsilon^2/8). Only |D| greater than this error
earns a nonzero result. At the center D=0, the relative
bound is undefined and stored as None, not infinity.
Nearby valid points and small masses can be inconclusive.
Crossing symmetry enforces the center cancellation of
the complete crossing-symmetric boundary difference;
a zero there does not test a whole analytic disc.

## Completeness at the first two-gauge onset

In this renormalizable candidate there is no tree
Phi^2 F^2 interaction or direct H Yukawa. A one-gauge
transition is color forbidden. A Phi single-particle
exchange is forbidden by the exact flavor-exchange
Phi parity. An H-to-two-gauge amplitude first needs
two loops: there is no direct H-fermion coupling.
Attaching a scalar contact or a heavy exchange to the
one-loop fermion box adds an internal scalar loop.
Mass, wavefunction and vertex counterterms also start
one order later. Thus the six finite-mass fermion boxes
are the full leading Phi Phi -> two-gauge transition.

Sewing two one-loop transitions with two gauge lines
gives 1+1+(2-1)=3 parent loops. A three-gauge final
state starts no earlier than 1+1+(3-1)=4 loops.
Scalar intermediate pairs start at invariant 4 and
heavy thresholds are outside the tested interval.
These facts identify the complete first two-gauge
discontinuity; they are not an enumeration of the
entire three-loop amplitude away from this channel.

Local counterterms have no discontinuity and cannot
remove the computed mismatch. The t-channel cut has
the same boundary prescription on both s boundaries.
Its forward limit is regular at this order: analyticity
of the massive transition and both Ward identities
supply one soft factor for each gauge leg, so the
two-gauge spectral term vanishes at least as t^2.
No Coulomb exchange is present for this neutral scalar.

The perturbative s/u cuts therefore still overlap the
old gapped disc even after this finite-mass correction.
This does not forbid choosing a local analytic germ
of logarithms, identify a nonperturbative confinement
spectrum or prove that every possible IR-subtracted
dispersion relation fails.
