# Exact state remainder and instantaneous projector comparison

Use E=sqrt(p^2+m0^2), m0=.99m, the same four exact frames
of S6.167, and real times. Their positive frequencies
satisfy e_j>=tau omega>=tau E. Put

    A=2pDelta/(tau E^3), r=1024/(tau E)<=1/1024.

The inherited complex-disk bounds imply
|theta_j|=|atan(g_j/e_j)|<=A r^j for j=0..3.
They control the full functions, not finite Taylor
approximations to the physical evolution.

In the fourth frame the exact state from either time
infinity differs from the instantaneous negative
projector by a transition amplitude no larger than

    integral_halfline |g4| <= C pDelta/(tau^4 E^6),
    C=2^42.

This is the exact unitary Duhamel inequality. Removing
the final diagonal phase leaves the negative projector
unchanged. The comparison projector is not asserted to
solve the Dirac equation or to define a Hadamard state.
The exact in/out states are the inherited Hadamard ones.

The four rotations in the physical mass eigenbasis
are Ux(theta0), Uy(theta1), Ux(theta2), Uy(theta3),
with the signs fixed in S6.167. Their Bloch z component is

    z4=cos0 cos1 cos2 cos3-cos0 sin1 sin3-sin0 sin2 cos3.

Thus its difference from cos0 is at most

    |theta0 theta2|+theta1^2+theta2^2+theta3^2
    <=A^2(2r^2+r^4+r^6)<3 A^2 r^2.

Here all cosines are positive and bounded by one,
|sin theta|<=|theta| and 1-cos theta<=theta^2/2.
There is no theta0 theta1 cross term: the alternating
axes matter. The mixed theta0 theta2 and theta1 theta3
terms are retained, not assumed zero.

One helicity in the first frame has physical energy

    -omega/sqrt(1+q0^2),
    q0=p Mdot/(2omega^3).

The physical energy operator in the mass eigenbasis is
omega sigma3, not the evolution generator with its
time-dependent basis connection g0 sigma2/tau. In
particular, neither -e1/tau nor -e4/tau replaces the
physical energy expectation. The connection remains
in the exact evolution and is not a stress term.

The absolute remainder after -omega+omega q0^2/2
is at most 3omega q0^4/8 by Taylor's integral formula
for (1+u)^(-1/2), u>=0. Its second derivative is at
most 3/4. Also |q0|<=pDelta/(2tau E^3).

Both occupied helicities therefore have the two UV terms

    -2omega + p^2 Mdot^2/(4omega^5).

The exact-state energy minus these terms is bounded
by the state-to-fourth, fourth-to-first and first-frame
Taylor remainders. No final connection or higher
transition order is discarded.
