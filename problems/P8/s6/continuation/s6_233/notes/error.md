# Uniform physical-angle tree error

The exact amplitude identity is proved in notes/amplitude.md. Let S=s>=4 and x=cos(theta) in[-1,1], with

t=-(S-4)(1-x)/2,
u=-(S-4)(1+x)/2.

Every shifted channel xi=channel-2 obeys abs(xi)<=S-2. For S<M_H²=D+2, every D-xi is positive. Therefore the exact remainder R is positive and

R<=3gamma(S-2)^4/[D-(S-2)].

The original tree is also positive on this domain. Its Galileon term3gamma S t u is nonnegative. Keeping only the s-channel lambda square and the full negative potential gives

A_original>=2lambda(S-2)²-8gamma.

The actual4lambda>8gamma and(S-2)²>=4 imply
A_original>=lambda(S-2)²>0.
The other lambda squares were retained as a positive discarded contribution, not neglected in the amplitude identity.

Since lambda=gamma D/2, division gives

0<R/A_original<=6r²/(1-r), r=(S-2)/D<1.

The right side is increasing because its derivative is6r(2-r)/(1-r)²>0. At S<=10^196, exact arithmetic gives

r<32/625,
6(32/625)²/(1-32/625)=6144/370625<1/60.

Thus the comparison is uniform over the entire physical angular interval and all massive energies2<=E<=10^98. It is not only a forward or central-angle estimate. The strict denominator and rational margins are checked independently.

This relative bound also controls the angular L2 error between these two TREE amplitudes. It does not identify the new tree with the full quantum amplitude appearing in S231. The parameter hierarchy and small tree values do not themselves control loop counterterms, physical mass/residue changes or resonance corrections.

The matching window is below the retained heavy pole. It is not a physical cutoff and does not remove the pole from a dispersive contour. A future physical matching calculation may have a different controlled window and must supply its own errors.

Numerical diagnostics at the actual hierarchy use1000 decimal digits: at low energies, contact/exchange pieces are hundreds of orders larger than their matched sum and its even smaller remainder. These diagnostics check the identities but are not the proof of the continuum bound.
