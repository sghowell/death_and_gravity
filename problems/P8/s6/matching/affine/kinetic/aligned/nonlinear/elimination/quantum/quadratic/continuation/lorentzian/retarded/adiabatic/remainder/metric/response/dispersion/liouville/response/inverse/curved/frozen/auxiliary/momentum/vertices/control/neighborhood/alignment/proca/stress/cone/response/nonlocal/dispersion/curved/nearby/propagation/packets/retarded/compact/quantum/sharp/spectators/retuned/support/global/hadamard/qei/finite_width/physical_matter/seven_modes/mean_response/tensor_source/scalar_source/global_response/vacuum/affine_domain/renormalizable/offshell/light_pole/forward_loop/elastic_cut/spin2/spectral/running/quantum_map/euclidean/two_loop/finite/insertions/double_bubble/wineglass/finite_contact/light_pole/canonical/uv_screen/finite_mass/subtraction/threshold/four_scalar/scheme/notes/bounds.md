# Exact logarithms, signed coefficient interval and complex pole bound

Use log(10)=3log(2)+log(5/4) with positive atanh
series ratios 1/3 and 1/9. The inherited log(2)
enclosure and a separate complete geometric tail
for log(5/4) give rational ell endpoints. At eight
retained terms they satisfy 900<ell_lower<ell_upper<1000.
The finite logarithm is specified, not replaced by
a guessed order-one coefficient.

Let T=4lambda, E_scalar the complete old one-loop
error, R=-2L+3g/D<0, E_box the S6.132 box bound,
fp_lower,fp_upper the S6.131 slope enclosure, and
0<r<r_upper the S6.112 scalar slope bound.
Using 144<Q<256, the extra relative term lies between

    R ell_upper/(2*144)+2(fp_lower-r_upper)-E_box/T
    R ell_lower/(2*256)+2fp_upper+E_box/T.

The negative sign of R reverses the correct log/Q
endpoints. These are not symmetric uncertainty
estimates silently dropping a finite scheme change.
Add the old interval [-E_scalar,E_scalar] to obtain
the complete one-loop coefficient interval.
The generic interface can return an inconclusive
nonpositive lower endpoint for valid wide errors.

At the actual parameters the extra relative effect
has magnitude below 10^-202. The inherited complete
scalar relative bound is approximately 4.470348358
times 10^-7; all new terms are retained as exact
rational enclosures even though this decimal display
cannot resolve them. Their sum is below 10^-6, and
the complete one-loop coefficient lower endpoint
is strictly positive.

For the pole, use the inherited scalar radius-two
coefficient g/[864 M^2(1-2/M)] and the unscaled
fermion bound B_F. Kappa is between
1-fp_upper and 1+r_upper-fp_lower. Therefore

    |f_F,R-Pi_scalar,R|/kappa
      <= (B_F+B_scalar)|s-1|^2/(1-fp_upper).

The coefficient is below 2 times 10^-405. Factoring
the inverse on |s-1|<=1 excludes any other zero,
while the exact anchored subtraction fixes the
mass and residue. The same estimate bounds the
negative curvature shift at s=0.

The scalar vacuum constant is negative, with absolute
value at most
(1+M^2)(ell_upper+3/2)/(4*144), below 10^396.
Subtract this upper bound from the fermion constant
lower endpoint, retaining both in the total reference.
The combined vacuum energy is between 10^799 and
10^800 before its fixed cancellation. No naturalness
prior or cosmological-constant exclusion is imposed.

All stated error budgets are for complete one-loop
quantities. They do not bound later primitive loops,
the exact gauge spectrum, high-energy contours or
finite-gravity/bounce matching.
