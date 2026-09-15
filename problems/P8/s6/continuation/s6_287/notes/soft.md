# Literal raw triangle sign and complete charge cancellation

The master convention is frozen by S283. The raw one-massless,
two-massive triangle in D=4+2EP is

C0mumu(t)=-Gamma(1-EP)*(4pi*nu^2)^(-EP)/(2EP)
          *integral_0^1[mu-t(1-v^2)/4]^(-1+EP)dv.

The minus sign comes from the odd number of Minkowski scalar
denominators; it is not reversed to make a Ward cancellation work.
Its raw IR pole is -2M(t)/EP. Expanding gives exactly the S283 finite
terms: -(1/EP+Gamma_E-log4pi)/2 times the inverse-parameter integral,
minus one-half its log(parameter/nu^2) integral.

Each scalar-graviton vertex is -i*T/sqrt(kappa). Three such vertices,
one graviton and two scalar propagators, and the loop integral
i*C0/(16pi^2), divided by the tree vertex -i*T/sqrt(kappa), give the
relative scalar-master coefficient -V_D/(16pi^2*kappa).
The literal soft numerator is

4[(p.pprime)^2-mu^2/(D-2)]
 =t^2-4mu*t+2mu^2+2mu^2*EP/(1+EP)
 =V_D.

Only the soft k0 numerator gives the IR pole; additional numerator
powers make that radial endpoint integrable. Thus the pole multiplying
the entire tree metric vertex is
+2V(t)M(t)/(16pi^2*kappa*EP).
This sign and normalization are independently checked against the
source vertex phases and the frozen triangle convention.

## All raw constants, not only a leading pole

The scalar inverse is K(p^2)=p^2-mu+Sigma(p^2), with S283's raw

Sigma_prime=mu/(16pi^2*kappa)
 [-3/EP+7+3log(4pi*nu^2/mu)-3Gamma_E].

Its UV part is -4mu/(16pi^2*kappa*EP), IR part
+mu/(16pi^2*kappa*EP), and the displayed finite remainder is retained.
In this pure-GR sector Sigma(mu)=0 in D. This is not a statement that
the whole theory has zero mass shift.

The background constant-metric identity gives
Gamma0=2pp*Kprime-eta*K.
Gauge.md proves the specified ordinary vertex equals it at the
on-shell point. Continuity.md proves that the ordinary F1 pole and
finite coefficients approach that value on the spacelike boundary.
Hence the proper F1 at zero transfer is1+Sigma_prime.

At a regulated formal simple pole, each scalar external leg contributes
Z^(1/2)=1-Sigma_prime/2 to first loop order. One endpoint has two such
halves, so its charge is
(1+Sigma_prime)(1-Sigma_prime)=1+O(hbar^2).
Two proper endpoints and all four scalar legs similarly cancel the
entire raw first-order correction. UV poles, IR poles, finite7,
Gamma_E and log4pi all cancel, not merely a selected logarithm.
Adding the same covariant kinetic UV/finite counterterm on both sides
preserves this equality; no new finite curvature coefficient is fixed.

## The finite-transfer remainder remains

After its two external-leg halves, a single endpoint still has pole
R(t)/(8pi^2*kappa*EP), where R=VM-mu/2.
The exact series is
R(t)=-11t/12+t^2/(10mu)+t^3/(84mu^2)+O(t^4).
For t=-tau the complete positive parameter identity is

R(-tau)=tau*integral_0^1
 [tau+mu(7+x^2)/2]/[4mu+tau(1-x^2)] dx.

For0<=tau<=mu it implies11tau/12<=R(-tau)<=7tau/6.
The lower comparison uses the denominator4mu+tau and its integrated
numerator; the upper uses4mu and then tau<=mu. The exact rational
margins are checked and agree with S286.

For fixed nonzero tau the pole does not vanish as EP tends to zero.
The result therefore gives continuity of its coefficient and the finite
coefficient, not a finite physical amplitude or permission to exchange
the transfer and infrared limits. Real/virtual inclusivity, dressing,
detector resolution and an exact interacting massless LSZ theorem
remain distinct questions.
