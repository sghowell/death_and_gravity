# Finite physical real-rate error

Fix the on-shell recoil map and pair-rest-frame angular bin in recoil.md.
The exact selected real integral is A0^-2 sum_helicities integral
d^3q/[(2pi)^3 2omega] J |M5|^2. The leading-soft integral uses J=1 and
A0 S0 evaluated at the exact Born on-shell momenta, with the same q.
Their common incoming flux and identical-final-scalar factor cancel.

Per helicity, write M5/A0=S0+delta, with |delta|<=B/sqrt(kappa),
B330000, |S0|<=64/(sqrt(kappa)omega),0<J<=1,1-J<=2omega.
The triangle inequality gives
|J|M5|^2/A0^2-|S0|^2| <=[B^2+(128B+8192)/omega]/kappa.
Integrate the entire sphere and both unit-normalized helicities:
the measure is omega d omega/(2pi^2). Thus for0<lambda<epsilon<=1/8,
the absolute difference is at most
[42248192epsilon+54450000000epsilon^2]/(2pi^2*kappa).
This dominates the difference integrably at zero, so its lambda0 limit
exists even though BOTH individual real rates still have a soft logarithm.

At epsilon1/8, pi>3 bounds the numerator divided by18 below10^8.
Original kappa10^800 gives strict error below10^-792, uniformly for
5/4<=E<=2,n>=128,mu1 and all angles. Multiplying by the positive
Born angular rate and integrating preserves the same relative bound.

Finite NumPy quadratures compare all graph tensors and both physical
polarizations independently: E1.25,n128,cos0.6; E1.5,n257,cos-0.3;
E2,n1024,cos0.2, at resolutions1/8,1/16, orders8/16. They calibrate
the implementation, not the proof or numerical integration error.
No numerical sign is extrapolated to a positivity claim.
