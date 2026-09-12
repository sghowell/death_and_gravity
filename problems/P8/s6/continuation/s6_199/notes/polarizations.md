# Complete physical stress amplitudes

Use eta=diag(+1,-1,-1,-1), canonical mass m>0, and E^2=k^2+m^2.
For momentum sign epsilon along e3, the creation readouts
(Electric,Magnetic,m A0,m Aspatial), with universal 1/sqrt(2E) omitted,
are

    T1 = (-i E e1, i epsilon k e2, 0, m e1)
    T2 = (-i E e2,-i epsilon k e1, 0, m e2)
    L  = (-i m epsilon e3,0,-k,epsilon E e3).

The omitted normalization is included exactly in the Lorentz-invariant
measure below. The longitudinal temporal component is constrained and
is not an additional oscillator.

For readouts u=(E,B,V0,V), v=(F,C,W0,W), the symmetric bilinear stress is

    T00 = (E.F+B.C+V0 W0+V.W)/2
    T0i = (E cross C+F cross B+V0 W+W0 V)_i/2
    Tij = -(E_i F_j+F_i E_j+B_i C_j+C_i B_j)/2
          +(V_i W_j+W_i V_j)/2
          +delta_ij (E.F+B.C+V0 W0-V.W)/2.

This is u^T M_ab v, not twice that amplitude. The connected Wick
exchange factor 2 is applied once, in the spectral measure.

Equivalently use the covariant Proca stress with
F_ab=-i(k_a A_b-k_b A_a), k.A=0 and A.A=-1. The tests construct it
directly from potentials, rather than from the displayed readout formula.
All nine pairs at +k and -k give symmetric T and T0mu=0 after exact
division by E^2-k^2-m^2. This proves the timelike Ward identity in the
center-of-mass frame; tensor covariance transfers it to every timelike
total momentum.

Completeness sum A_a A_b=-eta_ab+k_a k_b/m^2 is checked for both momenta.
Independent noncollinear rational Lorentz boosts retain the full Ward
identity and all nine amplitudes. No massless Maxwell substitution or
discarding of the longitudinal creation pair is involved.
