# Frozen physical inputs and attribution

S331 supplies positive Borel increments at fixed original hard E,u:
|deltaa_sigma-deltaa_tau|<=1700*T/kappa and
|deltaDelta_sigma-deltaDelta_tau|<=11000*T*(1+ln(1/T))/kappa
for sigma=tau+rho, mass(rho)=T and mass(sigma)<=1/8.
S325 gives the stronger Born bounds 1400*R/kappa and
10000*R*(1+ln(1/R))/kappa. The source module checks these actual APIs,
not duplicated fitted constants.

S330 proves the finite-cutoff renewal identity and relative conditioning
loss r=1-p/p_eta<=a*eta/x, retaining the zero-count atom eta^a.
S329 defines the fixed-Born leading probability and conditional moments,
using the earlier exact calorimetric CDF on energies <=1. We only use
that CDF inside the present x<=1/8 window. No extension above 1 is used.

S331 already proved the same-state regulator limit and cutoff estimates
for a and Delta separately. The new result controls the changed
remaining-energy logarithm and promotes the comparison to signed
densities on a specified common space. Its Delta density argument is
written anew: the earlier mean estimate alone is not an L1 proof.

The complete original trace, radial and amplitude-phase convention,
all primitive/matching records, six historical qualifications and
rejected S279 are retained. No unknown hard coefficient is assigned.
