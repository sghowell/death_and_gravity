# A finite spatial-frequency contribution of the original compact probes

The exact S6.91 multiplier has two leading fronts with amplitudes
below 10^6 and remainder bounded by C/k^2, C=10^44. Therefore,
for every k>=K0=4*10^31,

    |G_hat(t,s;k)| <= (2*10^6+C/K0)/k < B_G/k,
    B_G=10^13.

This is an absolute two-time bound for the same source and
observable, not a WKB relative estimate at a phase zero.
For the source high-spatial-frequency projection P_{k>Lambda},
Lambda>=K0, Plancherel gives

    |<f,G P_{k>Lambda} J>|
      <= B_G/Lambda ||f||_(L1 L2) ||J||_(L1 L2).

All physical-volume source factors have already been kept in
S6.91. The time integrations are over the original compact
supports, with retarded time order. Using their explicit norms
and rho=epsilon/20 gives

    ||f||_(L1 L2) ||J||_(L1 L2)
       <=64*20^(3/2)*a^2/epsilon
       <6400*a^2/epsilon.

Choose the exact finite comoving spatial cutoff

    Lambda=K0+25600*B_G*a^2/(D*epsilon).

It is above K0 and makes the high-spatial-frequency pairing
less than D/4. Since the full compact classical pairing is
at least 3D/4, its complementary contribution satisfies

    <f,G P_{k<=Lambda} J> >=D/2>0.

This is a frequency decomposition of the response to the
original compact probe. The projected source P_{k<=Lambda}J
is not compact; no simultaneous compactness and bandlimit
is claimed. It is the unprojected source and detector whose
full supports are matter-spacelike.

The cutoff is only spatial. Compact time bumps have unbounded
temporal frequency support, and no temporal tail/matching
bound has been supplied. Nor has Lambda been placed below
an interacting EFT or heavy-threshold cutoff. Restoring a
dimensionful time scale tau converts comoving momentum to
k/(tau e R); it does not, by itself, prove the missing
background, source, interaction or Wilson matching bounds.
