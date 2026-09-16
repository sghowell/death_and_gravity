# Positive Born amplitude, forward crossover and contour boundary

## Exact original Born forms

Use mass-one units and write Q=s-4,w=1-z^2. The complete tuned matter
Born amplitude is

A_m=g^2/(n-2)^2 sum_{a=s,t,u}(a-2)^2/(n-a),

where t=-Q(1-z)/2 and u=-Q(1+z)/2. This is positive for4<s<n.

The complete gravity Born amplitude, independently checked by full
four-vector stress contractions, is

kappa A_G=p/w+d+e w,
p=4Q+16+8/Q, d=-2s+6-2/s, e=Q^2/(4s).

For every Q>0,
p+2d=12+(4Q+32)/[Q(Q+4)]>0 and
d+e=-(7Q^2+40Q+40)/[4(Q+4)]<0. Thus

p/(2w)<kappa A_G<p/w

on0<w<=1. This proves positivity of the full Born denominator used in the
matter-loop interference bound, on its entire stated energy window.

## Explicit compact crossover

For25/4<=s<=16, n>=128,
257/9<=p<=194/3 and12<=sum(a-2)^2<=396.
All n-a lie between n-16 and n+12, and n-2>=63n/64.
Consequently4g^2/n^3<A_m<470g^2/n^3 and

xi/(33w)<A_G/A_m<17xi/w, xi=n^3/(kappa g^2).

Here all quantities are numerical in the chosen mass-one units.
At the exact original parameters,
1/(2*10^200)<xi<1/10^200.

For w>=10^-190, the ratio is below1.7*10^-9.
For w<=10^-204, it is above150. The existence of this inner cone does not
invalidate the selected matter soft calculation; it invalidates the
inference that large kappa makes every omitted Newton term uniformly small
at every angle. No angular cutoff is silently imposed as physical input.

## What a finite-contour Regge estimate would require

This is a conditional mathematical boundary statement, not a new
model-specific Regge bound.

At fixed negative t, put nu=s-2mu+t/2, so crossing is nu to-nu.
The original t-channel graviton numerator is
nu^2-2mu^2-t^2/4. Its nu^2 coefficient is-1/(kappa t).

Suppose an even analytic amplitude has been defined on the required
sheet, with an actual analytic neighborhood of nu0 after the explicitly
specified light s/u poles and cuts are accounted for. If further poles
remain, retain their residues explicitly. A finite-circle Cauchy deformation
then expresses the nu^2 coefficient as the signed cut integral plus those
pole contributions and

C_R(t)=(2pi i)^-1 integral_{|nu|=R} M(nu,t)dnu/nu^3.

The cut integral is NOT presumed positive at negative t. Positivity of its
forward limit needs the appropriate physical measure and limit control.

After the known s/u terms, the finite coefficient requires the combined
limit C_R(t)+1/(kappa t), not removal of the t pole in isolation.
For an actual Regge comparison M_R, the direct Cauchy estimate is

|C_R[M-M_R]| <= sup_{|nu|=R}|M-M_R|/R^2.

A relative approximation alone is insufficient near a1/t pole. For example,
on this finite contour only, let
M-M_R=nu^2 epsilon(1-exp(Lt))/(kappa t), epsilon,L>0.
The pole residue is unchanged and the normalized relative error on real
negative t is at most epsilon; its finite coefficient tends to
-epsilon L/kappa and has no L-independent bound.
This is a counterexample to that norm implication only, NOT a crossing-
complete, unitary or asymptotically admissible UV amplitude.

A genuinely sufficient alternative is a complex-transfer Cauchy bound:
if delta=kappa t(M-M_R)/nu^2 is holomorphic on and inside|t|=r,
delta(nu,0)=0 and its modulus on|t|=r is at most epsilon uniformly on
the energy contour, then |partial_t delta(nu,0)|<=epsilon/r.
The finite contour error is consequently at most epsilon/(kappa r).
Neither epsilon nor r is provided by the current low-energy graph bounds.

Moreover, original massless-pair cuts begin at s=0 at the next Newton
order. The complete amplitude does not inherit the tree's subthreshold
analytic disk merely from a massive external scalar. Its required
nonanalytic subtractions, observable dictionary and uniform limits must
be constructed, not inferred from a finite inclusive rate.

The primary FESR paper [Noumi and Tokuda](https://arxiv.org/abs/2212.08001)
provides a useful conditional route; its sections6.2-6.4 retain assumptions
on subleading terms and infrared null constraints, and its loop extension
does not establish this massive gravitational observable. Its numerical
bound is therefore not imported. The contour estimates above are the
explicit mathematical requirements used here, not a claimed application
of that paper to P8.
