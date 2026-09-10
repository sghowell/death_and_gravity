# Both one-loop sectors on a finite Euclidean reference window

Use the complete one-loop reference scheme of S6.133. Let r=Pi_scalar'(1)>0,
fp=f_fermion'(1)>0, and kappa=1+r-fp>0. The normalized selected inverse is

    Gamma(-t)=1+t+[f_fermion,R(-t)-Pi_scalar,R(-t)]/kappa.

For the old scalar kernel set Delta_1=x M+(1-x)^2, A=x(1-x), and
b_scalar=A/Delta_1. Its exact spacelike on-shell remainder is

    Pi_scalar,R(-t)=(g/Q) integral [d b_scalar-log(1+d b_scalar)] dx,
    d=t+1.

Thus 0<=Pi_scalar,R(-t)<=d r. Combine this with the fermion's negative
remainder, rather than dropping either sector or reversing its inverse
sign, to get

    0 <= 1-Gamma(-t)/(1+t)
       <= [r+C log(1+(t+1)/K)]/kappa,
    C=2 N Y/Q, K=4mF^2-1.

Now explicitly choose the test energy Lambda_ref=10^400 and
0<=t<=Lambda_ref^2. At mF=10^200, Lambda_ref=mF^2. Exact arithmetic gives

    mF^2-[1+(Lambda_ref^2+1)/K]=(3mF^4-5mF^2)/K>0.

The logarithm is consequently less than log(mF^2)=400 log(10).
Use the inherited positive-series rational enclosure of this specified
logarithm, the inherited scalar slope upper bound, the rational Yukawa
upper bound, 144<Q, and kappa>=1-fp_upper. The resulting defect is strictly
below 10^-202, so Gamma(-t)/(1+t) stays between 1-10^-202 and one.
In particular this selected one-loop inverse has no Euclidean zero on
that entire finite range. A wider valid error bound need not imply positivity.

The reference energy is a named mathematical test range. It is not a
Wilsonian cutoff inferred from this calculation, even if its numerical
value is also used as a Planck reference elsewhere. The flat-space kernel
here includes no gravitational corrections. Exact constant normalization
of the one-loop functional is not a later-loop error estimate or a claim
about reflection positivity or the complete Lorentzian spectrum.
