# Whole scalar cut and complete gravitational tree interference

## Full physical tree before angular integration

For s>4mu, put Q=s-4mu, t=-Q(1-z)/2,u=-Q(1+z)/2 and
D=4+2e. The exact minimal gravitational tree is
A_GR=-sum_a N_D(a;b,c)/(kappa*a),
N_D(a;b,c)=4mu^2(D-3)/(D-2)-2mu*a-b*c.
Combining all three exchanges, not only their endpoint poles, gives
kappa*A_GR=p/(1-z^2)+aa+bb*z^2,
p=4V_D(s)/Q,
aa=-7s/4+4mu+2mu^2/[s(1+e)],
bb=-Q^2/(4s).
The rational identity is checked before integration.

The normalized two-body angular measure has weight(1-z^2)^e.
Its exact moments are
average[1/(1-z^2)]=(1+2e)/(2e),
average[z^2]=1/(3+2e).
The singular numerator is integrated for positive e, before taking
any regulator limit. For independent numerical checks use the
convergent identity
integral_0^1(1-z^2)^(e-1)dz
=1/(2e)+integral_0^1(1-z^2)^e/(1+z)dz.
It follows by subtracting the endpoint in B(e,1/2)/2. It avoids
unresolved quadrature tails without discarding a finite term.

## Entire normal-sheet moment discontinuity

Only the s-channel massive moments have an imaginary part in the
physical region specified above. Their complete Feynman moments
are the light-root beta combination proved in S6.292, with its heavy
mass-square argument replaced by s. Equivalently the imaginary
part of integral(mu-s*x(1-x)-i0)^alpha is
-sin(pi*alpha)*(s/4)^alpha*beta^(2alpha+1)
 *B(1/2,alpha+1)/2,
beta=sqrt(1-4mu/s).
This is first used at alpha=e-1,e,e+1 for positive e; the full
analytic continuation and Coulomb phase are retained.

Gamma recurrence gives
ImJ/ImM=-2(1+2e)/(eQ),
ImM_(e+1)/ImM_e=-Q(e+1)/[2(2e+3)].
Since x(1-x)=(mu-A_s)/s, the entire weighted bubble ratio is
Yratio=mu/s+Q(e+1)/[2s(2e+3)].
Applying these identities to the whole-2T_e+E_e amplitude gives
the cut coefficient, after division by C*Phi2/(2kappa),
2V_D(s)(1+2e)/(eQ)-2(s-2mu)
 +(s+2mu/(1+e))*Yratio.
Its regular part equals aa+bb/(3+2e) exactly. Hence this is the
full gravitational tree average, including all finite regular terms.

## Phase and optical normalization

The complete dimension-dependent massive two-body phase is
Phi2=beta/(8pi)*(4pi nu^2)^(-e)*s^e*beta^(2e)
      *Gamma(1+e)/Gamma(2+2e).
The independent raw B0 discontinuity satisfies
ImB0/(16pi^2)=Phi2/2.
For identical intermediate scalars, the optical factor1/4 times
the two interference terms C*A_GR gives C*Phi2*average(A_GR)/2.
This agrees with the complete loop-moment discontinuity above.

The four independent60-digit numerical cut cases retain the complete
complex beta moments, the full radial Gamma/scale factors and the
convergently subtracted full angular tree. They test equality of
these whole expressions at positive e, not just agreement of poles.

## Scope of this independent check

At the isolated C/kappa coefficient there is no C tree transition
from two light scalars into two gravitons, nor an odd-matter
scalar-graviton transition. The selected one-loop nonlocal cut is
therefore the two-light-scalar interference above. Other original
coupling coefficients and cuts remain separate.

An interference coefficient is not a positive spectral density by
itself. A discontinuity cannot determine a local polynomial; the
tensor/on-shell calculation in notes/proper.md supplies that separate
part. This check proves neither full physical soft-factor unitarity
nor a detector/dressed observable, Regge control or original P8.
