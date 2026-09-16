# State-correct faces and the integrable remainder

The fixed-state leading soft theorem of S310 applies first at a=0
with b>0. Its hard state is sigma_b=(four massive recoiled legs,q2).
The same S300 map at a=0 is precisely the one-real recoil map.
The phase rho2 tends to rho1. Thus
 G(0,b)=b*S_a(sigma_b)*F1(b)/sqrt(kappa).
The other face is analogous. Taking the second leading limit removes
the marked-null contribution, which is O(b) in S_a, and both orders
give G00=S_a(Born)*S_b(Born)/kappa.

For a conserved signed state the leading current has contraction
with the soft momentum equal to the sum of the state's momenta.
The reduced massive legs alone sum to -q2; adding q2 gives zero.
The current is therefore Ward invariant only with the marked null
leg retained. Each projected term is also Ward invariant in the
other graviton through the inherited full47 amplitude. The Born
overlap and full434 amplitude retain their established Ward proofs.
This proves gauge consistency of this particular subtraction; it
does not assign separate physical meaning to individual graph classes.

The inclusion-exclusion definition subtracts both faces and restores
their common overlap exactly once:
 ab*R2=G-G(0,b)-G(a,0)+G00.
The complex estimate proves continuity of the faces at their common
origin. Apply the rectangle fundamental theorem first with positive
lower edges, then let them tend to zero. The integrable majorant
C/(u+v), C=10^-326, permits this limit and gives
 |R2|<C/(ab)*int_0^a du int_0^b dv/(u+v)=C*I(a,b)/(ab),
 I=(a+b)ln(a+b)-a ln a-b ln b.
If M=max(a,b),m=min(a,b),
 I=M ln(1+m/M)+m ln(1+M/m)
  <=m[1+ln(1+M/m)].
No assertion of a bounded pointwise hierarchical remainder is
needed; the possible logarithmic sharing envelope is integrable.

For the baseline, S311's complete current norm at a null total root
gives |F1(b)|<1/(64b), since sqrt(rho1)<=1. Each massive soft term
satisfies |pvec|^2/|p.n|<=12; four give48. For the marked null leg,
TT transversality bounds
 |q2.eps1.q2|/(q2.n1)<=b*sin(theta)^2/[1-cos(theta)]
 =b*(1+cos(theta))<=2b.
Thus a|J_a(sigma_b)|<49 uniformly in the relative angle, and
 |B2|<D/(ab),
 D=49/(32sqrt(kappa))+2304/kappa<2*10^-400.

Six exact original-parameter samples check both whole434 and47
faces, all massive shells, momentum conservation, correct-current
Ward identities and nonzero massive-only Ward failures.
Their amplitudes are rational before applying the exact phase roots.
Each unnormalized plus tensor has squared Frobenius norm2, so
every two-polarization expression is divided by2.
Rational triangle bounds, rho<=1 and I>min(a,b)/2 give exact
finite-point bound gates. Diagnostic numerical soft convergence is
not used as proof of the uniform estimate.
