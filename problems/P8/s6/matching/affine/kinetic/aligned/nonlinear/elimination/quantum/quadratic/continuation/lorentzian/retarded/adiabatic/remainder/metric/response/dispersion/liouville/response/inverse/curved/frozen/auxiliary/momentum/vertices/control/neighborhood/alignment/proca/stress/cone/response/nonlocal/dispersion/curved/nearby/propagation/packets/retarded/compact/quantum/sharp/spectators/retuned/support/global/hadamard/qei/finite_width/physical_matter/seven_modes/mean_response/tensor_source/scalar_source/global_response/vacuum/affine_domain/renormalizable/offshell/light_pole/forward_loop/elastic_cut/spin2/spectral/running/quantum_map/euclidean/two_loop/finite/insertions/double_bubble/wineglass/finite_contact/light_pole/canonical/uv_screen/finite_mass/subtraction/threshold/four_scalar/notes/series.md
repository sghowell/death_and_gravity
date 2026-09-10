# Every higher momentum degree is bounded

Use the frozen exact Euclidean propagator norm
||S0(q)||=1/sqrt(q^2+mF^2) and
||gamma.P||<=sum_mu |P_mu|. Expansion of each shifted
propagator is absolutely convergent in the domain of
the routing proof. For total insertion degree n,
the four propagators have binomial(n+3,3) weak
compositions. Some routes have a zero shift; counting
all four as variable is a safe overestimate.

After the entire local divergent reference has been
handled, every n>=1 term is absolutely integrable.
The common radial integral is

    integral d^4q/(2pi)^4 (q^2+mF^2)^(-2-n/2)
      =4/[Q n(n+2)mF^n].

The finite primitive and regulator boundary are inherited
from S6.129, whose immutable report is rebuilt.
Multiplying spin trace bound four, six cyclic orders,
six active color/flavor copies and radial factor four
gives the prefactor 576Y^2/Q. Scalar Yukawa vertices
have unit spin-matrix norm; there are no gauge
polarization factors in this calculation.

Put x=18/mF. For n>=4,

    binomial(n+3,3)/[n(n+2)]
      =(n+3)(n+1)/(6n) <= (n+5)/6.

Summing the majorant exactly gives

    sum_(n>=4) (n+5)x^n/6
      =x^4(9-8x)/[6(1-x)^2]
      <=(10/3)x^4 <4x^4  for 0<=x<=1/2.

The ratio derivative is (10-8x)/[6(1-x)^3]>0.
Small enumerated composition counts test the formula;
the generating-function argument bounds every term.
Odd degrees included in the estimate actually vanish,
which only makes it more conservative.

It follows that

    |R_F| <241864704Y^2/(QmF^4)
          <3 times 10^8 Y^2/(QmF^4).

Cauchy's coefficient estimate on |s-2|=1 bounds b2_F
by the same number. Lower momentum degrees are constant
on shell and contribute zero to this coefficient.
The estimate is not a claim about the sign of b2_F
alone or about omitted primitive loop orders.
