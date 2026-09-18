# Exact conditioned moments and retained original powers of kappa

Below reference energy1, the frozen leading CDF is
P(R<=r)=F(a)*r^a. For0<r<=x and a>0 its conditional CDF is
(r/x)^a, with density a*r^(a-1)/x^a. Direct integration gives
E[R^p|R<=x]=a*x^p/(a+p),p>0. At a=0 the energy is exactly
zero and all positive moments vanish, as does the continuous formula.

The retained S328 pointwise field bound is
||C_sigma-C_Born||sup<=23667R.
The Bochner triangle inequality and the first conditional moment imply
||E[C_sigma|cut]-C_Born||sup<=23667*a*x/(1+a).
Squaring the POINTWISE bound and using the second moment gives
E[||C_sigma-C_Born||sup^2|cut]<=23667^2*a*x^2/(2+a).
This is not the square of the mean. The full complex TT field is
kept, including its phase; no real-part replacement is made.

S299 proves a<4/(5kappa) on the original physical Born domain.
The exact positive roundups
23667*4/5<19000 and23667^2*2/5<230000000
give bounds19000*x/kappa and230000000*x^2/kappa before the
physical coefficient normalization. Dividing C by kappa^(3/2)
gives19000*x/kappa^(5/2) for the mean norm and
230000000*x^2/kappa^4 for the second moment.
At a=0 both quantities are exactly zero.

The extra smallness comes from the typical conditioned radiated energy,
not a changed physical parameter or a neglected phase. The unexpanded
calorimetric normalization cancels only in the exact conditional ratio.
These are selected marked leading-reference moments, not a full
radiative loop or a probability bound on the complete detector rate.
