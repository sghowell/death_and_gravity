# Restore the actual scalar mass with a fully convergent difference

The scalar boson in the actual model has mass one, not zero.
At fixed common MS parameters compare its full quadratic primitive
with the massless-boson reference. Route the chord as k-l, with no
external-soft momentum in that boson line. Its exact change is

 1/[(k-l)^2+1]-1/(k-l)^2
   =-1/{(k-l)^2[(k-l)^2+1]}.

The proper MS fermion kinetic/mass and Yukawa UV poles do not depend
on the boson mass. Their counterterm insertions into the first
fermion bubble are identical in the two calculations and cancel
in this difference. The whole four-fermion proper cycle gives only
a local mass tadpole after chord closure, hence has zero p^2 coefficient.
The overall p^2 UV pole is also mass independent.

Equivalently power counting the differentiated difference improves
the overall degree from zero to minus two. It improves the proper
fermion self-energy degree from one to minus one, and both proper
vertex degrees from zero to minus two. The whole-fermion-cycle
degree-zero local part is killed by the p^2 derivative. Thus the
difference is fully UV convergent, and the following joint bound
also checks its overlapping and infrared regions explicitly.

At each real k,l the three words have arc counts (3,1),(1,3),(2,2).
Put S_k=m^2+k^2,S_l=m^2+l^2. The coefficient of soft scaling degree
two is bounded by the circle supremum divided by R^2, with
R=min(sqrt(S_k),sqrt(S_l))/360. Each fermion resolvent is bounded
by 2/sqrt(S), so the trace majorant is 64/(S_k^(a/2)S_l^(b/2)).
This is one Taylor coefficient, not a full tail: no geometric factor
two is needed.

With x=(k-l)^2, use

 1/[x(x+1)] <=x^(-7/4),

because x^(3/4)<=1+x on both x<=1 and x>=1. Also
1/min(S_k,S_l)<=1/S_k+1/S_l. The resulting two terms per word
are exactly the convergent four-dimensional fractional sunset
integrals from S6.144, now with alpha+beta=3, gamma=7/4.
Their mass factor is (m^2)^(-3/4)=m^(-3/2).

All four numerator Gamma arguments lie between 1/4 and 9/4.
Euler's integral split at one bounds each by 4+6=10.
The three denominator Gamma values exceed 1/2 by the half-integer
recurrences. Each exact constant is therefore below 80000.
The three words, two minimum-radius terms and trace/coefficient
factor give

 |delta scalar MS slope|
 <=3981312000000 N Y_hi^2/(Q_lo^2 m^(3/2))
 < 4e12 N Y_hi^2/(Q_lo^2 m^(3/2)).

Every proper large-momentum region and the massless chord diagonal
are integrable. This majorant justifies taking the derivative,
integrating the difference and taking the regulator limit. There
is no missing finite epsilon times pole product in this already
convergent difference; those products were retained separately
in the massless MS reference. The bound does not assign a sign
to the finite mass correction.
