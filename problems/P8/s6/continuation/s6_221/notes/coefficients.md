# Actual fixed-profile coefficient bounds and principal pair

Use Theta=H-t/(1+t²)^4, E=1-3/[2(1+t²)^3], ell=1/[10(1+t²)^6], w=-ell E. Define
C=Theta E'-E Theta'+H E Theta-Theta² and F=C-w²/2.
The exact bare gradient is P(t²)/[800(1+t²)^18], where every nonzero coefficient of P is positive and P(0)=1199. The entire polynomial is in the report; no finite sample is used for its sign. Exactly,
J-F=1/[50(1+t²)^6].
At the bare bounce J=1215/800 and F/J=1199/1215.

For the fixed QG1 profile A=-pressure, B=-(rho+pressure)/2,
DeltaJ=(21delta²-3delta)A/2+(1-6delta)B and Tcorr=rho-3delta pressure.
Jc=J+DeltaJ. The original five-jet stress bound epsilon=1e-770 suffices for the zero-through-two jets used here. Product-rule envelopes in the certificate bound DeltaJ,Tcorr,A and their first two derivatives by1e-6; |DeltaJ|<4epsilon remains the sharper pointwise input.

On I, F>=1199/[800(5/4)^18]>1/100,
Jc>=1215/[800(5/4)^18]-4epsilon>1/100, and
Jc-F>=1/[50(5/4)^6]-4epsilon>1/1000.
The strict inequalities survive the unchanged retuning. Jc<100. The all-real bare identity is not extended to an unproved all-real retuned estimate.

For any rational profile with positive-even denominator, expand the numerator and bound its absolute coefficient sum with |t|<=1/2; bound the denominator below by its positive constant term. Applying this to each derivative gives:
Theta: (109/32,441/64,1057/64);
E: (93/64,9/2,99/4);
J: (3848766721573/107374182400,3993359652681/21474836480,241216786862451/214748364800);
ell: (1/10,3/5,51/10);
H: (2,5,13);
delta: (1/2,3/2,33/4).
These explicit real-axis envelopes and the product rule prove every coefficient jet through2 is below1e6. Sharper pointwise bounds used separately are |Theta|<=2, |H|<=2, |E|<=1, ell<=1/10, |w|<=1/10.

The full two-matter principal pair, after a valid chart congruence, is
Ktilde=[[2Jc+w²,w],[w,1]],
Gtilde=[[2F+w²,w],[w,1]].
Its characteristic determinant is2Jc(c²-1)(c²-F/Jc). Thus its speeds squared are1 andF/Jc, strictly separated on I. In the outer chart the congruence uses Theta; this formula alone is not a proof across Theta=0. The central chart's actual principal determinant is independently checked with E instead, so its characteristic statement remains valid through the crossing. No finite-q frozen-frequency argument is used.

For |t|<=1/4, E is negative, monotone in |t|, and lies between-1/2 and-1231/4913; hence1/4<|E|<=1/2. For |t|>=1/8, Theta=t[4-(1+t²)^-3]/(1+t²) gives |Theta|>=3/10. Both charts therefore exist at switches+-3/16. Since a²<=(5/4)^4, |P|>=100 implies q>=4096. The thresholds are analytical chart choices, not EFT cutoffs.
