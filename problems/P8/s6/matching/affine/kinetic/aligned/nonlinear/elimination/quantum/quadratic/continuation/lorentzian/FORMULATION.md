# P8-S6.63: physical-signature quadratic vector pole and evanescence

Date: 2026-09-08. Original P8 remains OPEN.

This checkpoint converts the local S6.61--62 quadratic counterterm
data to the prescribed Lorentzian matter metric. The action, state,
original clock history and lower-order matching remain unchanged.
It does not Wick rotate the selected state or replace the bounce
by a different time-rotated scale factor.

## Required local convention

Let Q_E denote the Euclidean loop residue in
Gamma_E,div=integral sqrt(g_E) Q_E/(32 pi^2 epsilon).
Local analytic continuation of the covariant jets gives

    H_E^(j)=(-i)^(j+1) H_L^(j),
    Y_E,t/s^(j)=(-i)^j Y_L,t/s^(j),
    k_E=k_L,
    Q_L=-Wick(Q_E).

The Yt,Ys amplitudes denote the temporal and spatial mixed mass
eigenvalue deviations. The physical covariant time-time component
is -Yt, not +Yt. The loop-pole sign is independently calibrated
against the frozen Lorentzian continuation of the Feynman integral.

The Lorentzian bare-action counterterm is -Q_L,D/(32 pi^2 epsilon).
Its normalized lapse variation is -E_L,3/(32 pi^2 epsilon)
plus the finite evanescent component 2(partial_D E_L,D)_3/(32 pi^2).
All D-dependent measure terms remain included.

## Acceptance gates

1. Construct the Lorentzian frame tensor derivatives and Einstein
   contractions independently, including the negative temporal metric.
2. Evaluate the auxiliary-metric heat coefficient with g00=-N^2
   directly and compare both derivative sectors with local analytic
   continuation as polynomial identities in arbitrary D.
3. Recover the frozen physical clock curvatures and all three
   analytically continued flat Feynman-pole coefficients.
4. Perform the actual lapse mass-profile projection and D-dependent
   compact variation in physical time; compare every coefficient
   and normalized operator through first dimensional order.
5. Keep an explicit nonzero control showing that naively replacing
   u by i*u in the scale factor changes the original history.
6. Pin all sources, fully rebuild S6.62, reject unsupported inputs
   and pass fresh ordinary and full P8 regressions.

## Scope

This supplies physical-signature local UV pole and specified
evanescent counterterm data. The existing full flat potential is
checked, not re-added. Every new derivative operator is quadratic
in Y and has zero value and first variation on the selected clock.

The Lorentzian finite bare/contact/retarded calculation, compatible
varied-state covariance, finite matching data and integrated
feedback remain open. No local pole sign supplies an exact
higher-derivative degree-of-freedom count, corrected quantum cone,
interacting cutoff or V/G/B verdict. Original P8 is not closed.
