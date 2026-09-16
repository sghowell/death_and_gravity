# P8 S6.299: complete leading-soft calorimetric sum

This continuation sums every multiplicity in the defined leading-soft
four-massive-scalar observable with a TOTAL unresolved energy cut. It
retains the complete S296 finite angular conversion and the original
analytic virtual convention. All original V/G/B/P8 gates remain OPEN.

The exact regulator-removed factor is
exp(Delta)*exp(-EulerGamma*a)/Gamma(1+a)*(E/nu)^a.
Its conversion remainder beyond1+Delta is below2*10^-1596 at the original
parameters, uniformly in0<E/nu<=1. The power of E/nu is never expanded
when its logarithm can scale with kappa.

See [formulation](FORMULATION.md), [whole sum](notes/poisson.md),
[conversion bound](notes/gamma.md), and [limit order](notes/limits.md).
This does NOT bound the full radiation amplitude minus its leading-soft
approximation uniformly in multiplicity. S295's single-real estimate
is not promoted to such a bound. Finite hard/Regge/bounce matching remains.

Native, direct, ordinary and CLI use original SymPy; the exact-GCD
adapter remains restricted to complete regression.
