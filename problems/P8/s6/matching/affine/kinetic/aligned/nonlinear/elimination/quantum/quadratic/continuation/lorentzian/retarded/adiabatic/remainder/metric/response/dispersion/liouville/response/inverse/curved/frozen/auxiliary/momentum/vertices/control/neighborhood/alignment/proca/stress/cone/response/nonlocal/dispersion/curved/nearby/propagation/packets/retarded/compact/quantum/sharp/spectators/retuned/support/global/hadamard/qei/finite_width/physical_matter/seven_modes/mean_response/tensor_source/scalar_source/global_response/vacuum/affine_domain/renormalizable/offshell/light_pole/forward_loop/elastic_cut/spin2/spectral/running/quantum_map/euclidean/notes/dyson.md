# Positive forms and only the known-kernel insertion series

The displayed one-loop Euclidean inverse is

    K(y)=(y+1)(1-Q(y)),  0<Q(y)<alpha<1.

Consequently (1-alpha)K0 <= K <= K0 as quadratic forms, where
K0=-Delta+1. These inequalities hold first for Schwartz test
fields and extend to the massive H1 form domain by the bounded
Fourier multiplier Q. The inverse satisfies

    G0 <= G <= G0/(1-alpha).

The corresponding source quadratic form is well-defined for
sources J with integral |Jhat(p)|^2/(1+p^2) dp finite. Its change
is at most alpha/(1-alpha) times the free source form. These
are Euclidean positivity statements, not Osterwalder-Schrader
reflection positivity or a reconstructed unitary quantum theory.

Expanding only the known multiplier gives sum_n Q^n. The norm
tail after retained order N is bounded by

    alpha^(N+1)/(1-alpha).

This is a convergent geometric family of repetitions of one
computed self-energy. It does not include primitive two-loop
self-energies, higher-loop counterterms or other diagram classes.

For two propagator lines, their possibly different ratios q,r
supply at total order n the sum over j=0,...,n of q^j r^(n-j).
There are n+1 distributions, each bounded by alpha^n. The tail
after N is

    alpha^(N+1)[N+2-(N+1)alpha]/(1-alpha)^2.

It follows by differentiating alpha^(N+2)/(1-alpha).
Native algebra checks the distinct-line coefficients and the
closed tails; pointwise tests cover the supported orders.
The software cap N<=32 is not a physical maximum loop order.

For an already absolutely integrable weighted two-line test
kernel, pointwise multiplication and domination give the same
relative tail bound against its absolute integral. An arbitrary
unsubtracted four-point trace is NOT such a test kernel.
For a local quartic insertion the ultraviolet radial behavior
includes

    integral_0^R y/(1+y)^2 dy
      = ln(1+R)+1/(1+R)-1,

which diverges logarithmically. Equivalently a generic nonzero
local potential sandwiched between two half-propagators is
marginal for the Hilbert-Schmidt norm in four dimensions.
The native primitive and asymptotic are a negative control
against applying a finite positive trace bound here.

Subtractions destroy the simple positive-integrand comparison
and require their own derivative and counterterm analysis.
Thus neither these tails nor the uniform propagator bound
establish a complete two-loop b2 error. They supply a rigorously
bounded building block for that separate task.
